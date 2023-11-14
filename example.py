
import mcmt
import open3d as o3d
import numpy as np
import time
import tqdm
import torch
import mcubes
import kaolin
from mesh_sdf import MeshSDF

num_sdf_query = 0

def sphere_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    return np.linalg.norm(points, axis=1) - 0.5

def empty_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    return -1


def box_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    return np.linalg.norm(points, axis=1) - 0.3

def lucy_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    mesh_name = "./assets/lucy_watertight.ply"
    sdf = True
    implicite_function = MeshSDF(mesh_name,true_sdf=sdf)
    return implicite_function(points)

def armadillo_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    # points = points - 0.5
    mesh_name = "armadillo.obj"
    sdf = True
    implicite_function = MeshSDF(mesh_name,true_sdf=sdf)
    return implicite_function(points)

def thai_statue_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    # points = points - 0.5
    mesh_name = "assets/thai_statue.obj"
    sdf = True
    implicite_function = MeshSDF(mesh_name,true_sdf=sdf)
    return implicite_function(points)


mesh_name = "./assets/vbunny.ply"

implicite_function = MeshSDF(mesh_name,true_sdf=True)


def happy_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]

    start_time = time.time()
    global implicite_function
    # print(f"MeshSDF Time: {time.time() - start_time}")
    return implicite_function(points)


if __name__ == "__main__":
    use_lloyd = True
    sdf_function = happy_sdf
    # sdf_function = thai_statue_sdf
    # sdf_function = empty_sdf 
    start_time = time.time()
    dim=1
    res = 256
    XX = np.linspace(-dim, dim, res)
    YY = np.linspace(-dim, dim, res)
    ZZ = np.linspace(-dim, dim, res)
    X, Y, Z = np.meshgrid(XX, YY, ZZ)

    points = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1) 
    sdfs = sdf_function(points)
    sdfs = sdfs.reshape(res, res, res).cpu().numpy()
    verts, faces = mcubes.marching_cubes(sdfs, 0)
    print(time.time() - start_time)
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(verts)
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    o3d.io.write_triangle_mesh("example_chair_gt_256.obj", mesh)
    exit()


    # sdf_function = sphere_sdf
    start_time = time.time()
    threshold = 1e-4
    resolution = 5
    counter = 0
    dim=2
    XX = np.linspace(-dim, dim, resolution)
    YY = np.linspace(-dim, dim, resolution)
    ZZ = np.linspace(-dim, dim, resolution)
    X, Y, Z = np.meshgrid(XX, YY, ZZ)
    points = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1)
    point_values = sdf_function(points)
    mcmt.add_points(points, point_values)
    # print(f"Total Query Number after initialization: {num_sdf_query}")

    # sampling random points
    for i in tqdm.tqdm(range(50)):
        sample_points = mcmt.sample_points_voronoi(128).reshape(-1, 3)
        # if use_lloyd:
        sample_points = mcmt.lloyd_relaxation(sample_points, 1, -0.5, 0.5).reshape(-1, 3)
        sample_values = sdf_function(sample_points)
        counter += sample_points.shape[0]
        mcmt.add_points(sample_points, sample_values)

    # for threshold in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1][::-1]:
    for i in tqdm.tqdm(range(100)):
        mid_points = mcmt.get_mid_points().reshape(-1, 3)
        mid_values = sdf_function(mid_points)
        next_points = mid_points[np.abs(mid_values) > threshold]
        # import pdb; pdb.set_trace()
        print(next_points.shape)
        counter += next_points.shape[0]
        if len(next_points) == 0:
            break
        mcmt.add_mid_points(next_points, mid_values[np.abs(mid_values) > threshold])
        # mcmt.output_triangle_mesh(f"example_chair_{i}.obj")

    print(f"MCMT Time: {time.time() - start_time}")

    print(counter)
    mcmt.output_triangle_mesh("example_chair.obj")
    # mcmt.output_grid_mesh("example_lucy_grid.obj", 0)


    # tet_verts = mcmt.get_grid_points().reshape(-1, 3)# - 0.5
    # tets = mcmt.get_grids().reshape(-1, 4)
    # sdf = sdf_function(tet_verts)
    # # import pdb; pdb.set_trace()
    # tet_verts = torch.tensor(tet_verts, dtype=torch.float32)
    # tets = torch.tensor(tets, dtype=torch.int64)
    # sdf = torch.tensor(sdf, dtype=torch.float32)
    # mesh_verts, mesh_faces = kaolin.ops.conversions.marching_tetrahedra(tet_verts.unsqueeze(0), tets, sdf.unsqueeze(0)) # running MT (batched) to extract surface mesh
    # mesh_verts = mesh_verts[0].cpu().numpy()
    # mesh_faces = mesh_faces[0].cpu().numpy()
    # mesh = o3d.geometry.TriangleMesh()
    # mesh.vertices = o3d.utility.Vector3dVector(mesh_verts)
    # mesh.triangles = o3d.utility.Vector3iVector(mesh_faces)
    # o3d.io.write_triangle_mesh("example123.obj", mesh)
    # print(f"Total Query Number after sampling: {num_sdf_query}")



    # counter = num_sdf_query
    # # optimization phase
    # start_time = time.time()
    # for i in tqdm.tqdm(range(100)):
    #     mid_points = mcmt.get_mid_points().reshape(-1, 3)
    #     mid_values = sdf_function(mid_points)
    #     # print(mid_values.shape)
    #     next_points = mid_points[np.abs(mid_values) > threshold]
    #     # print(next_points.shape)
    #     counter += next_points.shape[0]
    #     if len(next_points) == 0:
    #         break
    #     mcmt.add_mid_points(next_points, mid_values[np.abs(mid_values) > threshold])
    # print(time.time() - start_time)


    # # output
    # # print(f"Total Query Number after mid point optimization: {num_sdf_query}")
    # if use_lloyd:
    #     mcmt.output_triangle_mesh("example_lloyd.obj")
    # else:
    #     mcmt.output_triangle_mesh("example.obj")



    # points = mcmt.get_grid_points().reshape(-1, 3)
    # pcd = o3d.geometry.PointCloud()
    # pcd.points = o3d.utility.Vector3dVector(points)
    # o3d.visualization.draw_geometries([pcd])
