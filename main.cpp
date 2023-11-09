
//#include "CVT3D.h"
#include "fast_cvt3d.h"
#include <chrono>


void readPoints(std::string filename, std::vector<std::vector<int>> &sites) {
	std::ifstream infile(filename);
	float x, y, z;
	int quantized_x, quantized_y, quantized_z;
	std::string line;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		if (!(iss >> x >> y >> z)) { break; }
		quantized_x = (int)(x);  //[0 ~ resolution-1]
		quantized_y = (int)(y);  //[0 ~ resolution-1]
		quantized_z = (int)(z);  //[0 ~ resolution-1]
		std::vector<int> point = { quantized_x, quantized_y, quantized_z };
		sites.push_back(point);
	}
}

double sphere_sdf(double x, double y, double z)
{
	// Center of the sphere
	double centerX = 0.5;
	double centerY = 0.5;
	double centerZ = 0.5;

	// Radius of the sphere
	double radius = 0.4; // Adjust the radius as needed (should be less than 0.5 to fit within [0, 1])

	// Calculate the distance from the point to the center of the sphere
	double distance = std::sqrt((x - centerX) * (x - centerX) + (y - centerY) * (y - centerY) + (z - centerZ) * (z - centerZ));

	// Return the signed distance (negative inside the sphere, positive outside the sphere)
	return distance - radius;
}

int main(int argc, char **argv)
{
	using namespace GEO;

	int constrain_res = 4;
	int num_constrain_points = constrain_res * constrain_res * constrain_res;
	int num_points = 1000;
	int num_iter = 5;

	vector<double> points;
	points.resize(3 * (num_constrain_points + num_points));

	// always put constrained points in front
	for (int i = 0; i < constrain_res; i++)
	{
		for (int j = 0; j < constrain_res; j++)
		{
			for (int k = 0; k < constrain_res; k++)
			{
				double x = (1. / double(constrain_res - 1)) * i;
				double y = (1. / double(constrain_res - 1)) * j;
				double z = (1. / double(constrain_res - 1)) * k;
				int ind = i * constrain_res * constrain_res + j * constrain_res + k;
				points[ind * 3] = x;
				points[ind * 3 + 1] = y;
				points[ind * 3 + 2] = z;
			}
		}
	}

	for (index_t i = 0; i < 3 * num_points; ++i) {
		points[i + 3 * num_constrain_points] = Numeric::random_float64();
	}

	GEO::initialize();

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	Fast_CVT cvt = Fast_CVT();
	cvt.set_points(num_constrain_points + num_points, num_constrain_points, points.data());
	end = std::chrono::system_clock::now();
	std::cout << num_constrain_points + num_points << " points delaunay triangulation elapsed time: " << (end - start).count() / 1000 << " ms\n";

	cvt.output_delaunay("output.obj");

	std::vector<float> point_error;
	for (index_t i = 0; i < points.size() / 3; i++) {
		double sd = sphere_sdf(points[i * 3], points[i * 3 + 1], points[i * 3 + 2]);
		point_error.push_back(abs(1. / sd) * abs(1. / sd));
	}
	cvt.set_errors(point_error.size(), point_error.data());

	start = std::chrono::system_clock::now();
	cvt.Lloyd_iterations(num_iter);
	end = std::chrono::system_clock::now();
	std::cout << "CVT " << num_iter << " iters elapsed time: " << (end - start).count() / 1000 << " ms\n";

	cvt.output_delaunay("output2.obj");
	//int nb_points = 20000;
	//std::vector<vec3> samples = cvt.near_surface_sampling(nb_points);

	//std::ofstream outfile("sample_points.obj"); // create a file named "example.txt"
	//for (int i = 0; i < samples.size(); i++) {
	//	outfile << "v " << samples[i].x <<" "<< samples[i].y << " " << samples[i].z <<"\n";
	//}

	//std::vector<int> cell_indices;
	//for (int i = 0; i < 10000; i++) {
	//	cell_indices.push_back(1000);
	//}

	//std::vector<vec3> sample_points = cvt.sample_cells(cell_indices.data(), cell_indices.size());
	//std::ofstream outfile("sample_points.obj"); // create a file named "example.txt"

	//for (int i = 0; i < sample_points.size(); i++) {
	//	outfile << "v " << sample_points[i].x <<" "<< sample_points[i].y << " " << sample_points[i].z <<"\n";
	//}
	
}
