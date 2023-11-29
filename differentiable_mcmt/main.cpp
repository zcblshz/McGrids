#include "fast_mcmt.hpp"
#include "sdfs.hpp"

int main(int argc, char **argv)
{
    using namespace GEO;
    double threshold = 1e-5;
    int constrain_res = 32;
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

    points.clear();
    point_values.clear();
    for (int i = 0; i < constrain_res; i++)
    {
        for (int j = 0; j < constrain_res; j++)
        {
            for (int k = 0; k < constrain_res; k++)
            {
                double random_x = (double)rand() / RAND_MAX;
                double random_y = (double)rand() / RAND_MAX;
                double random_z = (double)rand() / RAND_MAX;
                points.push_back(Point(random_x-0.5, random_y-0.5, random_z-0.5));
                point_values.push_back(SDF::sphere_sdf(random_x-0.5, random_y-0.5, random_z-0.5));
            }
        }
    }
    // int num_vertices = mcmt.add_points(points, point_values);
    // std::cout << "Points added: " << num_vertices << std::endl;
    mcmt.export_grid_off("grid_parallel.off");
    mcmt.export_surface_obj("mesh.obj");

    // // // mcmt.clear();
    // std::vector<double> point_values;
    // for (int i = 0; i < points.size() / 3; i++)
    // {
    //     point_values.push_back(SDF::sphere_sdf(points[i * 3], points[i * 3 + 1], points[i * 3 + 2]));
    // }
    // mcmt.add_points(points.size() / 3, points.data(), point_values.data());
    // exit(0);
    // for(int i=0; i<0; i++){
    //     std::vector<double> sample_points = mcmt.sample_points_voronoi(256);
    //     point_values.clear();
    //     for (int i = 0; i < sample_points.size() / 3; i++)
    //     {
    //         point_values.push_back(SDF::sphere_sdf(sample_points[i * 3], sample_points[i * 3 + 1], sample_points[i * 3 + 2]));
    //     }
    //     mcmt.add_points(sample_points.size() / 3, sample_points.data(), point_values.data());
    //     std::cout << "Sampling iter: " << i << std::endl;
    //     mcmt.output_grid_points("grid_" + std::to_string(i) + ".obj");
    //     mcmt.save_triangle_mesh("mesh_" + std::to_string(i) + ".obj");


    // }


    // for (int i = 0; i < 00; i++)
    // {
    //     point_values.clear();

    //     std::string filename = "grid_" + std::to_string(i) + ".obj";
    //     std::ofstream mid_file(filename);
    //     // mcmt.save_triangle_mesh("before_mid_mesh_" + std::to_string(i) + ".obj");

    //     std::vector<double> mid_points = mcmt.get_mid_points();
    //     std::cout << "-------"  << std::endl;
    //     for (int j = 0; j < mid_points.size() / 3; j++)
    //     {
    //         point_values.push_back(SDF::sphere_sdf(mid_points[j * 3], mid_points[j * 3 + 1], mid_points[j * 3 + 2]));
    //         mid_file << "v " << mid_points[j * 3] << " " << mid_points[j * 3 + 1] << " " << mid_points[j * 3 + 2] << std::endl;
    //     }
    //     std::vector<double> additional_points;
    //     std::vector<double> additional_point_values;

    //     for(int j=0; j<point_values.size(); j++){
    //         if(std::abs(point_values[j]) > threshold){
    //             additional_points.push_back(mid_points[j*3]);
    //             additional_points.push_back(mid_points[j*3+1]);
    //             additional_points.push_back(mid_points[j*3+2]);
    //             additional_point_values.push_back(point_values[j]);
    //         }
    //     }
    //     if (additional_points.size() == 0)
    //         break;
    //     mcmt.add_mid_points(additional_points.size() / 3, additional_points.data(), additional_point_values.data());
    //     mcmt.save_triangle_mesh("step_2_mesh_" + std::to_string(i) + ".obj");
    //     mcmt.output_grid_points("step_2_grid_" + std::to_string(i) + ".obj");

    //     std::cout << "Mid Point iter: " << i << std::endl;
    // }

	// mcmt.clear();


    return 0;
}
