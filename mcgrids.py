import mcmt
import open3d as o3d
import numpy as np
import time
import tqdm
import torch
import mcubes
import kaolin
import os
import tempfile


class McGrids:

    def __init__(self, sdf_func, clip_min, clip_max, initial_resolution, num_sample_iters, num_sample_points, num_mid_iters, threshold, verbose=False, disable_cvt=False):
        self.sdf_func = sdf_func
        self.clip_min = clip_min
        self.clip_max = clip_max
        self.disbale_cvt = disable_cvt
        self.initial_resolution = initial_resolution
        self.num_sample_iters = num_sample_iters
        self.num_sample_points = num_sample_points
        self.num_mid_iters = num_mid_iters
        self.threshold = threshold
        self.verbose = verbose
        self.query_count = 0
        self.sdf_query_time = 0
        self.compute_time = 0
        self.completed = False
        mcmt.clear_mcmt()


    def __sdf__(self, points):
        self.query_count += points.shape[0]
        start_time = time.time()
        point_values = self.sdf_func(points)
        self.sdf_query_time += time.time() - start_time
        return point_values

    def __iters__(self):
        start_time = time.time()
        self.query_count = 0
        # create an a grid from uniform sampling
        XX = np.linspace(-self.clip_min, self.clip_max,
                         self.initial_resolution)
        YY = np.linspace(-self.clip_min, self.clip_max,
                         self.initial_resolution)
        ZZ = np.linspace(-self.clip_min, self.clip_max,
                         self.initial_resolution)
        X, Y, Z = np.meshgrid(XX, YY, ZZ)
        points = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1)
        point_values = self.__sdf__(points)
        mcmt.add_points(points, point_values)

        # sample from distribution and refine approximation
        for i in tqdm.tqdm(range(self.num_sample_iters), disable=not self.verbose):
            sample_points = mcmt.sample_points_voronoi(
                self.num_sample_points).reshape(-1, 3)
            if not self.disbale_cvt:
                sample_points = mcmt.lloyd_relaxation(sample_points, 1, -0.5, 0.5).reshape(-1, 3)
            sample_values = self.sdf_func(sample_points)
            mcmt.add_points(sample_points, sample_values)

        # mid point refinement
        pbar = tqdm.tqdm(range(self.num_mid_iters), disable=not self.verbose)
        for i in pbar:
            mid_points = mcmt.get_mid_points().reshape(-1, 3)
            mid_values = self.__sdf__(mid_points)
            pbar.set_description(f"MCMT: {mid_points.shape[0]}")
            next_points = np.ascontiguousarray(
                mid_points[np.abs(mid_values) > self.threshold])
            next_values = np.ascontiguousarray(
                mid_values[np.abs(mid_values) > self.threshold])
            if len(next_points) == 0:
                break
            if not self.disbale_cvt:
                next_points = mcmt.lloyd_relaxation(next_points, 1, -0.5, 0.5).reshape(-1, 3)
                next_values = self.__sdf__(next_points)
            mcmt.add_mid_points(next_points, next_values)
        self.completed = True
        self.compute_time = time.time() - start_time

    def extract_mesh(self):
        if not self.completed:
            self.__iters__()
            if self.verbose:
                print(f"MCMT Time: {self.compute_time}")
                print(f"SDF Query Count: {self.query_count}")
                print(f"Equalent to MarchingCube of resolution: {np.ceil(np.cbrt(self.query_count))} ")
                print(f"SDF Query Time: {self.sdf_query_time}")
        with tempfile.TemporaryDirectory() as tmpdirname:
            mcmt.output_triangle_mesh(os.path.join(tmpdirname, "mesh.obj"))
            mesh = o3d.io.read_triangle_mesh(
                os.path.join(tmpdirname, "mesh.obj"))
        vertices = np.asarray(mesh.vertices)
        faces = np.asarray(mesh.triangles)
        return vertices, faces

    def extract_grid(self):
        if not self.completed:
            self.__iters__()
        return mcmt.get_grid()
