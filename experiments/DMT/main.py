import torch
import kaolin
import tqdm
import time
import os
import numpy as np
from networks import VolumeSDF, Decoder
from scipy.spatial import Delaunay
from loss import loss_chamfer, laplace_regularizer_const, loss_eikonal
from utils import sample_points_from_obj_mesh
import argparse
import mcmt_helper

if __name__ == "__main__":
        
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--mesh_path', type=str, default='lucy_watertight.obj')
    argparser.add_argument('--network_type', type=str, default='decoder')
    argparser.add_argument('--assets_path', type=str, default='../../assets')
    argparser.add_argument('--device', type=str, default='cuda')
    argparser.add_argument('--lr', type=float, default=1e-4)
    argparser.add_argument('--laplacian_weight', type=float, default=0.1)
    argparser.add_argument('--iterations', type=int, default=10000)
    argparser.add_argument('--save_every', type=int, default=1)
    argparser.add_argument('--grid_res', type=int, default=128)
    argparser.add_argument('--num_sample_points', type=int, default=10000)
    argparser.add_argument('--use_dmt', action='store_true')
    argparser.add_argument('--additional_info', type=str, default='')

    args = argparser.parse_args()
    mesh_path = os.path.join(args.assets_path,  args.mesh_path)
    log_path = f"{args.mesh_path.split('.')[0]}_{args.network_type}_{'dmt' if args.use_dmt else 'mcmt'}_{args.additional_info}"
    timelapse = kaolin.visualize.Timelapse(log_path)
    points = sample_points_from_obj_mesh(mesh_path, args.num_sample_points).to(args.device)
    
    if points.shape[0] > 100000:
        idx = list(range(points.shape[0]))
        np.random.shuffle(idx)
        idx = torch.tensor(idx[:100000], device=points.device, dtype=torch.long)    
        points = points[idx]

    # The reconstructed object needs to be slightly smaller than the grid to get watertight surface after MT.
    points = kaolin.ops.pointcloud.center_points(points.unsqueeze(0), normalize=True).squeeze(0) * 0.9
    timelapse.add_pointcloud_batch(category='input', pointcloud_list=[points.cpu()], points_type = "usd_geom_points")

    # Initialize model and create optimizer
    if args.network_type == 'decoder':
        model = Decoder(multires=2).to(args.device)
        if os.path.exists(os.path.join(log_path, "model.pth")):
            model.load_state_dict(torch.load(os.path.join(log_path, "model.pth")))
        else:
            model.pre_train_sphere(1000)

    else:
        model = VolumeSDF().to(args.device)
        if os.path.exists(os.path.join(log_path, "model.pth")):
            model.load_state_dict(torch.load(os.path.join(log_path, "model.pth")))
        else:
            model.pre_train_sphere(1000)


    vars = [p for _, p in model.named_parameters()]
    optimizer = torch.optim.Adam(vars, lr=args.lr)
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lambda x: max(0.0, 10**(-x*0.0002))) # LR decay over time

    bar = tqdm.tqdm(range(args.iterations))
    tet_verts = None
    tets = None
    for it in bar:
        start_time = time.time()
        if it < 2000:
            if tet_verts is None:
                tet_verts_path = os.path.join(args.assets_path, '{}_verts.npz'.format(args.grid_res))
                tet_verts = torch.tensor(np.load(tet_verts_path)['data'], dtype=torch.float, device=args.device)
            if tets is None:
                tets = torch.tensor(([np.load(os.path.join(args.assets_path, '{}_tets_{}.npz'.format(args.grid_res, i)))['data'] for i in range(4)]), dtype=torch.long, device=args.device).permute(1,0)
        else:
            tet_verts, tets = mcmt_helper.get_tets(model, True, threshold=1e-4, max_iter=20)

        if args.network_type == 'decoder':
            pred = model(tet_verts) # predict SDF and per-vertex deformation
            sdf, deform = pred[:,0], pred[:,1:]
            mesh_verts, mesh_faces = kaolin.ops.conversions.marching_tetrahedra(tet_verts.unsqueeze(0), tets, sdf.unsqueeze(0)) # running MT (batched) to extract surface mesh
        else:
            sdf = model(tet_verts)
            mesh_verts, mesh_faces = kaolin.ops.conversions.marching_tetrahedra(tet_verts.unsqueeze(0), tets, sdf.unsqueeze(0)) # running MT (batched) to extract surface mesh
        mesh_verts, mesh_faces = mesh_verts[0], mesh_faces[0]
        # import open3d as o3d
        # mesh = o3d.geometry.TriangleMesh()
        # mesh.vertices = o3d.utility.Vector3dVector(mesh_verts.cpu().detach().numpy())
        # mesh.triangles = o3d.utility.Vector3iVector(mesh_faces.cpu().detach().numpy())
        # o3d.io.write_triangle_mesh("debug.ply", mesh)

        # pcd = o3d.geometry.PointCloud()
        # pcd.points = o3d.utility.Vector3dVector(points.cpu().detach().numpy())
        # o3d.io.write_point_cloud("debug_points.ply", pcd)

        # exit()
        chamfer_loss = loss_chamfer(mesh_verts, mesh_faces, points, 10000)
        laplace_loss = laplace_regularizer_const(mesh_verts, mesh_faces)
        loss = chamfer_loss + args.laplacian_weight * laplace_loss
        bar.set_description(f'Loss Total: {loss.item()} Loss Chamfer: {chamfer_loss.item()} Loss Laplacian: {laplace_loss.item()}  mesh vertices: {mesh_verts.shape[0]}, # of mesh faces: {mesh_faces.shape[0]}')
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()
        if (it) % args.save_every == 0 or it == (args.iterations - 1): 
            # save reconstructed mesh
            timelapse.add_mesh_batch(
                iteration=it+1,
                category='extracted_mesh',
                vertices_list=[mesh_verts.cpu()],
                faces_list=[torch.cat([mesh_faces.cpu(), torch.flip(mesh_faces.cpu(), dims=[-1])])+mesh_verts.cpu().max()]
            )
            # save parameters
            torch.save(model.state_dict(), os.path.join(log_path, "model.pth"))
