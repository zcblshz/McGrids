#include "fast_mcmt.hpp"
#include "sdfs.hpp"

int main(int argc, char **argv)
{
    double threshold = 1e-6;
    int constrain_res = 8;
    int num_constrain_points = constrain_res * constrain_res * constrain_res;

    std::vector<Point> points;
    std::vector<double> point_values;
    std::vector<double> point_curvatures;


    for (int i = 0; i < constrain_res; i++)
    {
        for (int j = 0; j < constrain_res; j++)
        {
            for (int k = 0; k < constrain_res; k++)
            {   
                double px = (double)i / (constrain_res - 1);
                double py = (double)j / (constrain_res - 1);
                double pz = (double)k / (constrain_res - 1);

                double random_x = (double)rand() / RAND_MAX * 1e-6;
                double random_y = (double)rand() / RAND_MAX * 1e-6;
                double random_z = (double)rand() / RAND_MAX * 1e-6;
                px = px + random_x;
                py = py + random_y;
                pz = pz + random_z;
                double value = SDF::sphere_sdf(px-0.5, py-0.5, pz-0.5);
                points.push_back(Point(px-0.5, py-0.5, pz-0.5));
                point_values.push_back(value);
                point_curvatures.push_back(1);

            }
        }
    }



    MCMT::MCMT mcmt = MCMT::MCMT();
    int num_initial_points = mcmt.add_points(points, point_values, point_curvatures);

    // mcmt.export_grid_off("Grid_initial.off");
    // mcmt.export_surface_obj("Surface_initial.obj");

    std::cout << "Initial points added: " << num_initial_points << std::endl;

    // sample iters
    for(int i =0; i< 100; i++){
        std::vector<Point> sampled_points = mcmt.sample_tetrahedron(128);
        std::vector<double> sampled_point_values;
        std::vector<double> mid_point_curvatures;
        for (int j = 0; j < sampled_points.size(); j++)
        {
            double value = SDF::sphere_sdf(sampled_points[j].x(), sampled_points[j].y(), sampled_points[j].z());
            sampled_point_values.push_back(value);
            mid_point_curvatures.push_back(1);
        }
        int num_vertices = mcmt.add_points(sampled_points, sampled_point_values, mid_point_curvatures);
        std::cout << "Iter: " << i << " Total Points: " << num_vertices << std::endl;
    }

    // mcmt.export_grid_off("Grid_after_random_sample.off");
    // mcmt.export_surface_obj("Surface_after_random_sample.obj");

    // mid point iters
    for(int i=0; i< 100; i++){
        std::cout <<    "Iteration: " << i << std::endl;
        std::vector<Point> mid_points = mcmt.get_mid_points();
        std::cout <<    "Num mid points: " << mid_points.size() << std::endl;
        std::vector<Point> added_points;
        std::vector<double> mid_point_values;
        std::vector<double> mid_point_curvatures;

        for (int i = 0; i < mid_points.size(); i++)
        {
            double value = SDF::sphere_sdf(mid_points[i].x(), mid_points[i].y(), mid_points[i].z());
            if(std::abs(value) > threshold){
                added_points.push_back(mid_points[i]);
                mid_point_values.push_back(value);
                mid_point_curvatures.push_back(1);
            }
        }
        if (added_points.size() == 0)
            break;
        int num_vertices = mcmt.add_points(added_points, mid_point_values, mid_point_curvatures);
        std::cout << "Total Points: " << num_vertices << std::endl;
    }

    mcmt.export_grid_off("Grid_after_mid_point.off");
    mcmt.export_surface_obj("Surface_after_mid_point.obj");

    return 0;
}
