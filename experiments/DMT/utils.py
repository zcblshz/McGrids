import kaolin
import numpy as np
import torch

def read_obj_mesh(mesh_path):
    mesh = kaolin.io.obj.import_mesh(mesh_path)
    verts = mesh.vertices
    faces = mesh.faces
    verts = verts - verts.min(dim=0)[0]
    verts = verts / verts.max()
    verts -= 0.5
    return verts, faces

def sample_points_from_obj_mesh(mesh_path, num_points):
    verts, faces = read_obj_mesh(mesh_path)
    points = kaolin.ops.mesh.sample_points(verts.unsqueeze(0), faces, num_points)[0][0]
    return points