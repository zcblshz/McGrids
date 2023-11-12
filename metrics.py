import torch
import torch.nn.functional as F
import kaolin
import numpy as np
from experiments.DMT.utils import *
import open3d as o3d
from scipy.spatial import cKDTree
import os

import torch
from scipy.spatial import cKDTree

def compute_normal_difference(gt_vertices, gt_faces, input_vertices, input_faces, num_points=100000):

    # Sample points from meshes
    input_points, input_indices = kaolin.ops.mesh.sample_points(input_vertices.unsqueeze(0), input_faces, num_points)
    #
    # Step 1. Calculate normals for sample points on input mesh
    #
    # Convert tensors to numpy arrays and flatten to 2D
    input_points_np = input_points.cpu().numpy().reshape(-1, 3)
    # calculate face normals
    input_face_normals = torch.cross(input_vertices[input_faces[:, 1]] - input_vertices[input_faces[:, 0]],
                                    input_vertices[input_faces[:, 2]] - input_vertices[input_faces[:, 0]], dim=1)
    input_face_normals = F.normalize(input_face_normals, p=2, dim=1)
    sample_points_normal_input = input_face_normals[input_indices]

    #
    # Step 2. Calculate normals for sample points on gt mesh
    #
    # Build KD trees
    gt_kd_tree = cKDTree(gt_vertices)
    # Query KD trees to find nearest neighbors
    _, gt_indices = gt_kd_tree.query(input_points_np, k=3)

    # Convert indices to PyTorch tensors
    gt_indices = torch.from_numpy(gt_indices).to(input_points.device)
    # Calculate face normals
    gt_face_normals = torch.cross(gt_vertices[gt_faces[:, 1]] - gt_vertices[gt_faces[:, 0]],
                                    gt_vertices[gt_faces[:, 2]] - gt_vertices[gt_faces[:, 0]], dim=1)
    # Calculate the vertex normals
    vertex_normals = torch.zeros((gt_vertices.shape[0], 3), dtype=gt_vertices.dtype, device=gt_vertices.device)
    # Expand the face_normals tensor to match the shape of face_tensor for element-wise addition
    expanded_normals = gt_face_normals[gt_faces.view(-1)].view(*gt_faces.shape, -1)
    # Use index_add_ to accumulate the normals for each vertex
    vertex_normals.index_add_(0, gt_faces.view(-1), expanded_normals.view(-1, 3))
    # Normalize the accumulated normals to get average normals
    vertex_normals = F.normalize(vertex_normals, p=2, dim=1)

    sample_points_normal_gt = vertex_normals[gt_indices].sum(1)
    # Normalize the accumulated normals to get average normals
    sample_points_normal_gt = F.normalize(sample_points_normal_gt, p=2, dim=1)

    # Calculate normal difference
    normal_difference = 1 - torch.cosine_similarity(sample_points_normal_gt, sample_points_normal_input, dim=2)

    # Average the differences globally
    average_difference = normal_difference.mean()

    return average_difference

def compute_edge_chamfer(gt_vertices, gt_faces, input_vertices, input_faces, num_points=100000):

    # Sample points from input mesh
    input_points, input_indices = kaolin.ops.mesh.sample_points(input_vertices.unsqueeze(0), input_faces, num_points)
    input_points_np = input_points.cpu().numpy().reshape(-1, 3)
    # calculate face normals
    input_face_normals = torch.cross(input_vertices[input_faces[:, 1]] - input_vertices[input_faces[:, 0]],
                                    input_vertices[input_faces[:, 2]] - input_vertices[input_faces[:, 0]], dim=1)
    input_face_normals = F.normalize(input_face_normals, p=2, dim=1)
    input_sample_normals = input_face_normals[input_indices].squeeze()

    input_kd_tree = cKDTree(input_points_np)
    # Query KD trees to find nearest neighbors
    neighbor_indices = input_kd_tree.query_ball_point(input_points_np, 0.01, p=2.)

    output = []
    for neighbor_index in neighbor_indices:
        neighbor_index = torch.from_numpy(np.array(neighbor_index)).to(input_face_normals.device).type(torch.long)
        input_neighbor_normals = input_sample_normals[neighbor_index]
        expanded_self_normal = input_neighbor_normals[0,:].expand(input_neighbor_normals.size(0), -1)
        normal_similarity = torch.cosine_similarity(expanded_self_normal, input_neighbor_normals, dim=1)
        # Check if the minimum normal similarity is smaller than 0.1
        if torch.min(normal_similarity) < 0.1:
            output.append(neighbor_index[0].item())

    output = torch.from_numpy(np.asarray(output)).type(torch.long)
    # p = input_points[torch.from_numpy(output)]
    # Convert the output list to a tensor

    input_sample_tensor = torch.tensor(input_points[0][output]).to('cuda')

    # Sample points from gt mesh
    gt_points, gt_indices = kaolin.ops.mesh.sample_points(gt_vertices.unsqueeze(0), gt_faces, num_points)
    gt_points_np = gt_points.cpu().numpy().reshape(-1, 3)
    # calculate face normals
    gt_face_normals = torch.cross(  gt_vertices[gt_faces[:, 1]] - gt_vertices[gt_faces[:, 0]],
                                    gt_vertices[gt_faces[:, 2]] - gt_vertices[gt_faces[:, 0]], dim=1)
    gt_face_normals = F.normalize(gt_face_normals, p=2, dim=1)
    gt_sample_normals = gt_face_normals[gt_indices].squeeze()

    gt_kd_tree = cKDTree(gt_points_np)
    # Query KD trees to find nearest neighbors
    neighbor_indices = gt_kd_tree.query_ball_point(gt_points_np, 0.01, p=2.)

    output = []
    for neighbor_index in neighbor_indices:
        neighbor_index = torch.from_numpy(np.array(neighbor_index)).to(gt_face_normals.device).type(torch.long)
        gt_neighbor_normals = gt_sample_normals[neighbor_index]
        expanded_self_normal = gt_neighbor_normals[0, :].expand(gt_neighbor_normals.size(0), -1)
        normal_similarity = torch.cosine_similarity(expanded_self_normal, gt_neighbor_normals, dim=1)
        # Check if the minimum normal similarity is smaller than 0.1
        if torch.min(normal_similarity) < 0.1:
            output.append(neighbor_index[0].item())
    output = torch.from_numpy(np.asarray(output)).type(torch.long)

    # Convert the output list to a tensor
    gt_sample_tensor = torch.tensor(gt_points[0][output]).to('cuda')
    chamfer_dist = kaolin.metrics.pointcloud.chamfer_distance(input_sample_tensor.unsqueeze(0), gt_sample_tensor.unsqueeze(0))

    #mesh = o3d.geometry.TriangleMesh()
    #mesh.vertices = o3d.utility.Vector3dVector(gt_sample_tensor[0].detach().cpu())
    #o3d.io.write_triangle_mesh("data/gt_points_edge.obj", mesh)
    #mesh.vertices = o3d.utility.Vector3dVector(gt_sample_tensor[0].detach().cpu())
    #o3d.io.write_triangle_mesh("data/input_points_edge.obj", mesh)

    return chamfer_dist


if __name__ == "__main__":

    gt_mesh_path = "data/armadillo_marching_cube_128.obj"
    input_mesh_path = "data/armadillo_marching_cube_64.obj"

    # Read in the meshes using kaolin
    gt_vertices, gt_faces = read_obj_mesh(gt_mesh_path)
    input_vertices, input_faces = read_obj_mesh(input_mesh_path)

    print("evaluating the models.")

    # sample on mesh
    num_points = 100000
    gt_points = kaolin.ops.mesh.sample_points(gt_vertices.unsqueeze(0), gt_faces, num_points)[0]
    input_points = kaolin.ops.mesh.sample_points(input_vertices.unsqueeze(0), input_faces, num_points)[0]
    gt_points_cuda = gt_points.to('cuda')
    input_points_cuda = input_points.to('cuda')
    chamfer_dist = kaolin.metrics.pointcloud.chamfer_distance(gt_points_cuda, input_points_cuda)
    print("Chamfer Distance:", chamfer_dist.item())

    #mesh = o3d.geometry.TriangleMesh()
    #mesh.vertices = o3d.utility.Vector3dVector(gt_points[0].cpu())
    #o3d.io.write_triangle_mesh("data/gt_points.obj", mesh)
    #mesh.vertices = o3d.utility.Vector3dVector(input_points[0].cpu())
    #o3d.io.write_triangle_mesh("data/input_points.obj", mesh)

    normal_consistency = compute_normal_difference(gt_vertices, gt_faces, input_vertices, input_faces, num_points=10000)

    print("Normal Consistency:", normal_consistency.item())

    edge_chamfer_dist = compute_edge_chamfer(gt_vertices, gt_faces, input_vertices, input_faces, num_points=100000)
    print("Edge Chamfer Distance:", edge_chamfer_dist.item())


