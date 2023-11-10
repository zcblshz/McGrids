
import mcmt
import open3d as o3d
import numpy as np
import time
import tqdm
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
    use_lloyd = False
    sdf_function = lucy_sdf
    threshold = 1e-4
    resolution = 5
    XX = np.linspace(0, 1, resolution)
    YY = np.linspace(0, 1, resolution)
    ZZ = np.linspace(0, 1, resolution)
    X, Y, Z = np.meshgrid(XX, YY, ZZ)
    points = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1)
    point_values = sdf_function(points)
    mcmt.add_points(points, point_values)
    # print(f"Total Query Number after initialization: {num_sdf_query}")

    # sampling random points
    for i in tqdm.tqdm(range(20)):
        sample_points = mcmt.sample_points(256).reshape(-1, 3)
        if use_lloyd:
            sample_points = mcmt.lloyd_relaxation(sample_points, 1).reshape(-1, 3)
        sample_values = sdf_function(sample_points)
        mcmt.add_points(sample_points, sample_values)
    # print(f"Total Query Number after sampling: {num_sdf_query}")

    counter = num_sdf_query
    # optimization phase
    start_time = time.time()
    for i in tqdm.tqdm(range(100)):
        mid_points = mcmt.get_mid_points().reshape(-1, 3)
        mid_values = sdf_function(mid_points)
        # print(mid_values.shape)
        next_points = mid_points[np.abs(mid_values) > threshold]
        # print(next_points.shape)
        counter += next_points.shape[0]
        if len(next_points) == 0:
            break
        mcmt.add_points(next_points, mid_values[np.abs(mid_values) > threshold])
    print(time.time() - start_time)


    # output
    # print(f"Total Query Number after mid point optimization: {num_sdf_query}")
    if use_lloyd:
        mcmt.output_triangle_mesh("example_lloyd.obj")
    else:
        mcmt.output_triangle_mesh("example.obj")



    points = mcmt.get_grid_points().reshape(-1, 3)
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points)
    o3d.visualization.draw_geometries([pcd])
