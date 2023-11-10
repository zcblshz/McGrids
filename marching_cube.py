
import mcmt
import open3d as o3d
import numpy as np
import time
import tqdm
import mcubes
from mesh_sdf import MeshSDF

num_sdf_query = 0

def sphere_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    points = points - 0.5
    return np.linalg.norm(points, axis=1) - 0.3

def box_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    points = points - 0.5
    return np.linalg.norm(points, axis=1) - 0.3

def lucy_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    points = points - 0.5
    mesh_name = "lucy_watertight.ply"
    sdf = True
    implicite_function = MeshSDF(mesh_name,true_sdf=sdf)
    return implicite_function(points)

def armadillo_sdf(points):
    global num_sdf_query
    num_sdf_query += points.shape[0]
    points = points - 0.5
    mesh_name = "armadillo.obj"
    sdf = True
    implicite_function = MeshSDF(mesh_name,true_sdf=sdf)
    return implicite_function(points)


if __name__ == "__main__":
    sdf_function = lucy_sdf
    resolution = 512
    start_time = time.time()
    XX = np.linspace(0, 1, resolution)
    YY = np.linspace(0, 1, resolution)
    ZZ = np.linspace(0, 1, resolution)
    X, Y, Z = np.meshgrid(XX, YY, ZZ)
    points = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1)
    point_values = -sdf_function(points)
    point_values = np.asarray(point_values)
    vertices, triangles = mcubes.marching_cubes(point_values.reshape(resolution, resolution, resolution), 0)
    vertices = np.concatenate([vertices[:, 1:2], vertices[:, :1], vertices[:, 2:]], axis=1)
    triangles =  np.concatenate([triangles[:, 2:], triangles[:, 1:2], triangles[:, :1]], axis=1) # [N,3]

    vertices = ((vertices) / (resolution-1))

    print(f"Total Query Number: {points.shape[0]}")
    print(f"Took {time.time() - start_time} seconds")
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(vertices)
    mesh.triangles = o3d.utility.Vector3iVector(triangles)
    o3d.io.write_triangle_mesh(f"lucy_marching_cube_{resolution}.ply", mesh)

