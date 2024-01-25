import os
import mcubes
import torch
import torch.nn as nn
import open3d as o3d
import numpy as np
# import fastsweep
# import drjit
import tqdm
# from diff_operators import *

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class SDF(nn.Module):

    def __init__(self):
        super(SDF, self).__init__()

    def forward(self, points):
        raise NotImplementedError

    def normals(self, points, batch_size=2 ** 16):
        # points [N, 3]
        points.requires_grad = True
        # sdf = self.sdf(points)
        normals = []
        for batch_points in torch.split(points, batch_size):
            sdf = self.sdf(batch_points)
            normals.append(torch.autograd.grad(sdf, batch_points, grad_outputs=torch.ones_like(
                sdf), create_graph=True, retain_graph=True, only_inputs=True)[0].cpu().detach())
        normal = torch.cat(normals, dim=0)
        # return torch.nn.functional.normalize(normal, dim=-1)
        return normal

    def sdf(self, points):
        raise NotImplementedError

    def __grid(self, resolution, size):
        X = torch.linspace(-size / 2., size / 2., resolution)
        Y = torch.linspace(-size / 2., size / 2., resolution)
        Z = torch.linspace(-size / 2., size / 2., resolution)
        XX, YY, ZZ = torch.meshgrid(X, Y, Z, indexing='ij')
        points = torch.stack([XX, YY, ZZ], dim=-1).reshape(-1, 3)
        return points.to(device)

    def __marching_cubes(self, sdfs, resolution, size, iso, redistance=False):
        sdfs = sdfs.reshape(resolution, resolution, resolution).cpu().detach().numpy()
        if redistance:
            data = drjit.cuda.TensorXf(sdfs)
            sdfs = fastsweep.redistance(data).numpy()
        vertices, triangles = mcubes.marching_cubes(sdfs, iso)
        vertices = (((vertices + 0.5) / (resolution)) - 0.5)  # [-0.5 ~ 0.5]
        vertices *= size / ((resolution - 1) / resolution)
        # flip triangle winding
        triangles = triangles[:, [2, 1, 0]]
        return vertices, triangles

    def export_mesh(self, file_name, resolution=16, size=2, iso=0):
        vertices, triangles = self.compute_mesh(resolution, size, iso)
        mesh = o3d.geometry.TriangleMesh()
        mesh.vertices = o3d.utility.Vector3dVector(vertices)
        mesh.triangles = o3d.utility.Vector3iVector(triangles)
        #vertex_normals = self.normals(torch.from_numpy(vertices).float().to(device)).cpu().detach().numpy()
        #mesh.vertex_normals = o3d.utility.Vector3dVector(vertex_normals)
        # create output folder if not exists
        os.makedirs(os.path.dirname(os.path.abspath(file_name)), exist_ok=True)
        o3d.io.write_triangle_mesh(file_name, mesh)

    def compute_mesh(self, resolution, size=2, iso=0, batch_size=2 ** 16):
        points = self.__grid(resolution=resolution, size=size)
        sdfs = []
        for batch_points in torch.split(points, batch_size):
            sdfs.append(self.sdf(batch_points).cpu().detach())
        sdfs = torch.cat(sdfs)
        vertices, triangles = self.__marching_cubes(sdfs, resolution, size, iso)
        return vertices, triangles.astype(np.int32)

    def compute_mesh_redistance(self, resolution, size=2, iso=0, batch_size=2 ** 16):
        points = self.__grid(resolution=resolution, size=size)
        sdfs = []
        norms = []
        for batch_points in torch.split(points, batch_size):
            sdfs.append(self.sdf(batch_points).cpu().detach())
            norms.append(self.normals(batch_points).norm(dim=-1).cpu().detach())
        sdfs = torch.cat(sdfs)
        norms = torch.cat(norms)
        sdfs = sdfs / norms
        vertices, triangles = self.__marching_cubes(sdfs, resolution, size, iso)
        return vertices, triangles.astype(np.int32)

    def compute_sdf_and_curvature(self, points):
        # points [N, 3]
        points = torch.from_numpy(points).float().to(device)
        points.requires_grad = True
        sdf = self.sdf(points)
        return sdf.cpu().detach().numpy(), mean_curvature(sdf, points).cpu().detach().numpy()

class ExtrudeSDF(SDF):

    def __init__(self, height=0.3):
        super(ExtrudeSDF, self).__init__()
        height = torch.tensor(height)
        self.height = nn.Parameter(height, requires_grad=False)

    def forward(self, points):
        return self.sdf(points)

    def sdf(self, points):
        # points [N, 3]
        # height [1]
        extrude_height = torch.abs(self.height)
        sdf_2d = self.sdf_2d(points[..., :2])
        h = torch.abs(points[..., -1]) - extrude_height.unsqueeze(0)
        d = torch.stack([sdf_2d, h], dim=-1)
        return d[..., 0].max(d[..., 1]).clamp_max(
            0) + torch.norm(d.clamp_min(0), dim=-1)

    def sdf_2d(self, points):
        raise NotImplementedError


class ExtrudeCylinder(ExtrudeSDF):

    def __init__(self, radius=0.3, height=0.3):
        super(ExtrudeCylinder, self).__init__(height=height)
        self.radius = radius

    def sdf_2d(self, points):
        # points [N, 2]
        return torch.norm(points, dim=-1) - self.radius


class ExtrudeHexagram(ExtrudeSDF):

    def __init__(self, radius=0.3, height=0.3):
        super(ExtrudeHexagram, self).__init__(height=height)
        radius = torch.tensor(radius)
        self.radius = nn.Parameter(radius, requires_grad=False)
        k = torch.tensor([-0.5, 0.8660254038, 0.5773502692, 1.7320508076])
        self.k = nn.Parameter(k, requires_grad=False)

    def sdf_2d(self, points):
        # points [N, 2]
        p = torch.abs(points)
        p = p - 2.0 * torch.min(torch.einsum("d,id->i",
                                             self.k[:2],
                                             p),
                                torch.zeros_like(p[:,
                                                 0])).unsqueeze(-1) * self.k[:2].unsqueeze(0)
        p = p - 2.0 * torch.min(torch.einsum("d,id->i", self.k[[1, 0]], p), torch.zeros_like(
            p[:, 0])).unsqueeze(-1) * self.k[[1, 0]].unsqueeze(0)
        p = p - torch.stack([torch.clamp(p[...,
        0],
                                         self.radius * self.k[2],
                                         self.radius * self.k[-1:]),
                             self.radius.repeat(points.shape[0])],
                            dim=-1)
        return torch.norm(p, dim=-1) * torch.sign(p[..., 1])


class SphereSDF(SDF):

    def __init__(self, radius=0.3):
        super(SphereSDF, self).__init__()
        self.radius = radius

    def forward(self, points):
        return self.sdf(points)

    def sdf(self, points):
        # points [N, 3]
        return torch.norm(points, dim=-1) - self.radius

    def curvature(self, points):
        sdf = self.sdf(points)
        mean_curvature = self.mean_curvature(sdf, points)
        return mean_curvature

    def mean_curvature(self, y, x):
        grad = self.gradient(y, x)
        grad_norm = torch.norm(grad, dim=-1)
        unit_grad = grad.squeeze(-1)/(grad_norm.unsqueeze(-1))

        Km = 0.5*self.divergence(unit_grad, x)
        return Km

    def gradient(self, y, x, grad_outputs=None):
        if grad_outputs is None:
            grad_outputs = torch.ones_like(y)
        grad = torch.autograd.grad(
            y, [x], grad_outputs=grad_outputs, create_graph=True
        )[0]
        return grad
    
    def divergence(self, y, x):
        div = 0.
        for i in range(y.shape[-1]):
            div += torch.autograd.grad(
                y[..., i], x, torch.ones_like(y[..., i]), create_graph=True
            )[0][..., i:i+1]
        return div

class BoxSDF(SDF):

    def __init__(self, dims=torch.tensor([0.5, 0.5, 0.5])):
        super(BoxSDF, self).__init__()
        self.dims = dims
        self.dims = nn.Parameter(self.dims, requires_grad=False)

    def forward(self, points):
        return self.sdf(points)

    def sdf(self, points):
        # points [N, 3]
        # dims [3]
        N, _ = points.shape
        dims = torch.abs(self.dims)
        q_points = points.abs() - dims.unsqueeze(0).repeat(N, 1)
        lengths = (q_points.max(torch.zeros_like(q_points))).norm(dim=-1)
        zeros_points = torch.zeros_like(lengths)
        xs = q_points[..., 0]
        ys = q_points[..., 1]
        zs = q_points[..., 2]
        filling = ys.max(zs).max(xs).min(zeros_points)
        return lengths + filling


class LinkSDF(SDF):

    def __init__(self, le=0.3, r1=0.3, r2=0.2):
        super(LinkSDF, self).__init__()
        le = torch.tensor(le)
        r1 = torch.tensor(r1)
        r2 = torch.tensor(r2)
        self.le = nn.Parameter(le, requires_grad=False)
        self.r1 = nn.Parameter(r1, requires_grad=False)
        self.r2 = nn.Parameter(r2, requires_grad=False)

    def forward(self, points):
        return self.sdf(points)

    def sdf(self, points):
        q = torch.stack([points[:, 0], torch.clamp(torch.abs(points[:, 1]) - self.le, min=0.0), points[:, 2]], dim=-1)
        return (torch.stack([(q[:, :2]).norm(dim=-1) - self.r1, q[:, 2]], dim=-1)).norm(dim=-1) - self.r2


class RoundedBoxSDF(SDF):

    def __init__(self, dims=torch.tensor([0.2, 0.2, 0.2]), r=0.1):
        super(RoundedBoxSDF, self).__init__()
        self.dims = dims
        self.r = torch.tensor(r)
        self.dims = nn.Parameter(self.dims, requires_grad=False)
        self.r = nn.Parameter(self.r, requires_grad=False)

    def forward(self, points):
        return self.sdf(points)

    def sdf(self, points):
        q = torch.abs(points) - self.dims.unsqueeze(0)
        return torch.clamp(q, min=0.0).norm(dim=-1) + torch.clamp(
            torch.maximum(q[:, 0], torch.maximum(q[:, 1], q[:, 2])), max=0.0) - self.r

    def curvature(self, points):
        sdf = self.sdf(points)
        mean_curvature = self.mean_curvature(sdf, points)
        return mean_curvature

    def mean_curvature(self, y, x):
        grad = self.gradient(y, x)
        grad_norm = torch.norm(grad, dim=-1)
        unit_grad = grad.squeeze(-1)/(grad_norm.unsqueeze(-1)+1e-4)

        Km = 0.5*self.divergence(unit_grad, x)
        return Km

    def gradient(self, y, x, grad_outputs=None):
        if grad_outputs is None:
            grad_outputs = torch.ones_like(y)
        grad = torch.autograd.grad(
            y, [x], grad_outputs=grad_outputs, create_graph=True
        )[0]
        return grad
    
    def divergence(self, y, x):
        div = 0.
        for i in range(y.shape[-1]):
            div += torch.autograd.grad(
                y[..., i], x, torch.ones_like(y[..., i]), create_graph=True
            )[0][..., i:i+1]
        return div
    
if __name__ == "__main__":
    SDF = RoundedBoxSDF()
    points = torch.randn(1000, 3)
    points.requires_grad = True
    curvature = SDF.curvature(points)[:, 0]
    import plyfile
    # write ply with curvature
    vertex = np.array(points.detach().cpu())
    vertex_color = np.array(curvature.detach().cpu())
    vertex_normal = np.array(SDF.normals(points).detach().cpu())
    vertex = np.hstack((vertex, vertex_color[:, None]))
    vertex = np.hstack((vertex, vertex_normal))
    vertex = np.core.records.fromarrays(vertex.transpose(), names='x, y, z, curv, nx, ny, nz')
    el = plyfile.PlyElement.describe(vertex, 'vertex')
    plydata = plyfile.PlyData([el], text=True)
    plydata.write('curvature.ply')
    