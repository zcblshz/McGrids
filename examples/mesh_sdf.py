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
        unsigned_distance = self.scene.compute_distance(query_point)
        signed_distance = self.scene.compute_signed_distance(query_point)
        return signed_distance.numpy()

    def __call__(self, points):
        return self.sdf(points)