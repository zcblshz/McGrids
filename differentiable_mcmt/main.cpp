#include "fast_mcmt.hpp"
#include "sdfs.hpp"

int main(int argc, char **argv)
{
    using namespace GEO;
    double threshold = 1e-5;
    int constrain_res = 8;
    int num_constrain_points = constrain_res * constrain_res * constrain_res;

    std::vector<Point> points;
    std::vector<double> point_values;

    // always put constrained points in front
    for (int i = 0; i < constrain_res; i++)
    {
        for (int j = 0; j < constrain_res; j++)
        {
            for (int k = 0; k < constrain_res; k++)
            {
                double random_x = (double)rand() / RAND_MAX;
                double random_y = (double)rand() / RAND_MAX;
                double random_z = (double)rand() / RAND_MAX;
                double value = SDF::sphere_sdf(random_x-0.5, random_y-0.5, random_z-0.5);
                points.push_back(Point(random_x-0.5, random_y-0.5, random_z-0.5));
                point_values.push_back(value);
            }
        }
    }

    MCMT mcmt = MCMT();
    int num_initial_points = mcmt.add_points(points, point_values);
    std::cout << "Initial points added: " << num_initial_points << std::endl;


    for(int i=0; i< 20; i++){
        std::cout <<    "Iteration: " << i << std::endl;
        std::vector<Point> mid_points = mcmt.get_mid_points();
        std::cout <<    "Num mid points: " << mid_points.size() << std::endl;
        std::vector<Point> added_points;
        std::vector<double> mid_point_values;
        for (int i = 0; i < mid_points.size(); i++)
        {
            double value = SDF::sphere_sdf(mid_points[i].x(), mid_points[i].y(), mid_points[i].z());
            if(std::abs(value) > threshold){
                added_points.push_back(mid_points[i]);
                mid_point_values.push_back(value);
            }
        }
        if (added_points.size() == 0)
            break;
        int num_vertices = mcmt.add_points(added_points, mid_point_values);
        std::cout << "Total Points: " << num_vertices << std::endl;
        // mcmt.export_grid_off("Grid_" + std::to_string(i) + ".off");
        // mcmt.export_surface_obj("Surface" + std::to_string(i) + ".obj");
    }

    mcmt.export_grid_off("Grid.off");
    mcmt.export_surface_obj("Surface.obj");

    return 0;
}
