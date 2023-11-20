import argparse
import sys
import os
import open3d as o3d
from nglod_sdf import NGLODSDF
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lib.models import *
from lib.options import parse_options
from mcgrids import McGrids



if __name__ == "__main__":


    parser = parse_options(return_parser=True)
    app_group = parser.add_argument_group('app')
    app_group.add_argument('--output_path', type=str,
                        required=True, help='Path to output mesh')
    app_group.add_argument('--threshold', type=float,
                        default=1e-4, help='Terminating threshold')
    app_group.add_argument('--resolution', type=int,
                        default=4, help='Initial resolution')
    app_group.add_argument('--num_sample_iters', type=int,
                        default=20, help='Number of sample iterations')
    app_group.add_argument('--num_sample_points', type=int,
                        default=128, help='Number of sample points')
    app_group.add_argument('--num_mid_iters', type=int,
                        default=300, help='Number of mid iterations')
    app_group.add_argument('--min', type=float, default=1, help='min clip')
    app_group.add_argument('--max', type=float, default=1, help='max clip')
    args = parser.parse_args()

    implicite_function = NGLODSDF(args)
    grids = McGrids(implicite_function, clip_min=args.min, clip_max=args.max, initial_resolution=args.resolution, num_sample_iters=args.num_sample_iters,
                            num_sample_points=args.num_sample_points, num_mid_iters=args.num_mid_iters, threshold=args.threshold, verbose=True)
    vertices, faces = grids.extract_mesh()
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(vertices)
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    o3d.io.write_triangle_mesh(args.output_path, mesh)
