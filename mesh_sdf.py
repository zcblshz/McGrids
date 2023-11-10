import torch
import numpy as np
import open3d as o3d
import matplotlib.pyplot as plt

class MeshSDF:

    def __init__(self, mesh_file_path, true_sdf):
        self.mesh_file_path = mesh_file_path
        self.mesh = o3d.io.read_triangle_mesh(self.mesh_file_path)
        self.mesh = o3d.t.geometry.TriangleMesh.from_legacy(self.mesh)
        # normalize mesh to unit cube
        vertices = self.mesh.vertex.positions.cpu().numpy()
        center = vertices.min(axis=0)
        scale = np.max(np.max(vertices, axis=0) - np.min(vertices, axis=0))
        self.mesh.translate(-center)
        self.mesh.scale(scale=1/scale, center=[0,0,0])
        self.mesh.translate(-self.mesh.get_center())
        self.mesh.scale(scale=0.8, center=[0, 0, 0])
        self.true_sdf = true_sdf
        self.scene = o3d.t.geometry.RaycastingScene()
        _ = self.scene.add_triangles(self.mesh) 

    def sdf(self, points):
        if len(points.shape) == 1:
            points = points.reshape(1, -1)
        query_point = o3d.core.Tensor(points, dtype=o3d.core.Dtype.Float32)

        # Compute distance of the query point from the surface
        unsigned_distance = self.scene.compute_distance(query_point)
        signed_distance = self.scene.compute_signed_distance(query_point)
        return signed_distance
    

    def sdf_2d_grid(self, resolution):
        XX = np.linspace(-1, 1, resolution)
        YY = np.linspace(-1, 1, resolution)
        X, Y = np.meshgrid(XX, YY)
        points = np.hstack((X.reshape(-1, 1), Y.reshape(-1, 1)))
        points = np.hstack((points[:, :1], np.zeros((points.shape[0], 1)), points[:, 1:]))
        sdf = self.sdf(points).cpu().numpy()
        sdf = sdf.reshape(resolution, resolution)
        # plot zero level set        
        plt.imshow(sdf)
        self.mesh.translate(-self.mesh.get_center())
        plt.show()
        return sdf

    def __call__(self, points):
        if points.shape[-1] == 2:
            points = np.hstack((points[:, :1], np.zeros((points.shape[0], 1)), -points[:, 1:]))
        sdf = self.sdf(points).cpu().numpy()
        if self.true_sdf:
            return torch.from_numpy(sdf)
        else:
            return torch.sigmoid(torch.from_numpy(sdf*-10)) - 0.5


if __name__ == "__main__":
    mesh_sdf = MeshSDF("tiktok.obj")
    points = np.random.randn(1000,2)
    sdf = mesh_sdf.sdf_2d_grid(resolution=1024)
    import pdb; pdb.set_trace()
    # print(mesh_sdf(torch.from_numpy(points)))