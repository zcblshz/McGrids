import argparse
import sys
import os
import open3d as o3d
from mesh_sdf import MeshSDF
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcgrids import McGrids



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', type=str,
                        required=True, help='Path to input mesh')
    parser.add_argument('--output_path', type=str,
                        required=True, help='Path to output mesh')
    parser.add_argument('--threshold', type=float,
                        default=1e-3, help='Terminating threshold')
    parser.add_argument('--resolution', type=int,
                        default=4, help='Initial resolution')
    parser.add_argument('--num_sample_iters', type=int,
                        default=20, help='Number of sample iterations')
    parser.add_argument('--num_sample_points', type=int,
                        default=128, help='Number of sample points')
    parser.add_argument('--num_mid_iters', type=int,
                        default=100, help='Number of mid iterations')
    parser.add_argument('--min', type=float, default=0.5, help='min clip')
    parser.add_argument('--max', type=float, default=0.5, help='max clip')
    args = parser.parse_args()

    implicite_function = MeshSDF(args.input_path, true_sdf=True)
    grids = McGrids(implicite_function, clip_min=args.min, clip_max=args.max, initial_resolution=args.resolution, num_sample_iters=args.num_sample_iters,
                            num_sample_points=args.num_sample_points, num_mid_iters=args.num_mid_iters, threshold=args.threshold, verbose=True)
    vertices, faces = grids.extract_mesh()
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(vertices)
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    o3d.io.write_triangle_mesh(args.output_path, mesh)
