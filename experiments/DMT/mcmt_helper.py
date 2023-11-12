from importlib import reload
import mcmt
import tqdm
import torch
import numpy as np
from scipy.spatial import Delaunay

def sdf_function(points, model, device="cuda"):
    points = points - 0.5
    sdf =  model(torch.tensor(points, dtype=torch.float32, device=device)).cpu().detach().numpy()
    if len(sdf.shape) == 2:
        sdf = sdf[:, 0]
    return sdf


resolution = 8
XX = np.linspace(0, 1, resolution)
YY = np.linspace(0, 1, resolution)
ZZ = np.linspace(0, 1, resolution)
X, Y, Z = np.meshgrid(XX, YY, ZZ)
points = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1)

def get_tets(model, optimize_mid_points, threshold=1e-3, max_iter=10, device="cuda"):
    mcmt.clear_mcmt()
    use_lloyd = True
    point_values = sdf_function(points, model, device=device)
    mcmt.add_points(points, point_values)
    for i in range(10):
        sample_points = mcmt.sample_points_rejection(512).reshape(-1, 3)
        if use_lloyd:
            sample_points = mcmt.lloyd_relaxation(sample_points, 1).reshape(-1, 3)
        sample_values = sdf_function(sample_points, model, device=device)
        mcmt.add_points(sample_points, sample_values)
    if optimize_mid_points:
        for i in tqdm.tqdm(range(max_iter)):
            mid_points = mcmt.get_mid_points().reshape(-1, 3)
            mid_values = sdf_function(mid_points, model, device=device)
            next_points = mid_points[np.abs(mid_values) > threshold]
            if len(next_points) == 0:
                break
            mcmt.add_mid_points(next_points, mid_values[np.abs(mid_values) > threshold])
    tet_verts = mcmt.get_grid_points().reshape(-1, 3) - 0.5
    tets = mcmt.get_grids().reshape(-1, 4)
    return torch.tensor(tet_verts, dtype=torch.float32, device=device), torch.tensor(tets, dtype=torch.long, device=device)