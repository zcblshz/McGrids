#include "mcmt.hpp"
#include <tbb/tbb.h>

namespace GEO
{
	MCMT::MCMT()
	{
		delaunay_ = new PeriodicDelaunay3d(periodic_, 1.0);
		if (!periodic_)
		{
			delaunay_->set_keeps_infinite(true);
		}
	}

	MCMT::~MCMT()
	{
		delete delaunay_;
	}

	void MCMT::add_points(int num_points, double *point_positions, double *point_values)
	{
		for (index_t i = 0; i < num_points; i++)
		{
			point_positions_.push_back(point_positions[i * 3]);
			point_positions_.push_back(point_positions[i * 3 + 1]);
			point_positions_.push_back(point_positions[i * 3 + 2]);
			point_values_.push_back(point_values[i]);
		}
		delete delaunay_;
		delaunay_ = new PeriodicDelaunay3d(periodic_, 1.0);
		if (!periodic_)
		{
			delaunay_->set_keeps_infinite(true);
		}
		delaunay_->set_vertices(point_positions_.size() / 3, point_positions_.data());
		delaunay_->compute();
	}

	std::vector<double> MCMT::interpolate(double *point1, double *point2, double sd1, double sd2)
	{
		double p1_x =  point1[0];
		double p1_y =  point1[1];
		double p1_z =  point1[2];
		double p2_x =  point2[0];
		double p2_y =  point2[1];
		double p2_z =  point2[2];
		double t = sd1 / ((sd1 - sd2));
		return std::vector<double>{p1_x + t * (p2_x - p1_x), p1_y + t * (p2_y - p1_y), p1_z + t * (p2_z - p1_z)};
	}

	std::vector<double> MCMT::compute_face_mid_point(int num_points, const std::vector<double> &points)
	{
		double mid_point_x = 0;
		double mid_point_y = 0;
		double mid_point_z = 0;

		for (int i = 0; i < num_points; i++)
		{
			mid_point_x += points[i * 3];
			mid_point_y += points[i * 3 + 1];
			mid_point_z += points[i * 3 + 2];
		}
		mid_point_x /= num_points;
		mid_point_y /= num_points;
		mid_point_z /= num_points;
		return std::vector<double>{mid_point_x, mid_point_y, mid_point_z};
	}



	void MCMT::save_face(std::ofstream &output_mesh, const std::vector<double> &points, int &vertex_count)
	{
		int point_size = points.size() / 3;

		if (point_size == 3)
		{
			output_mesh << "v " << points[0*3] << " " << points[0*3+1] << " " << points[0*3+2] << "\n";
			output_mesh << "v " << points[1*3] << " " << points[1*3+1] << " " << points[1*3+2] << "\n";
			output_mesh << "v " << points[2*3] << " " << points[2*3+1] << " " << points[2*3+2] << "\n";
			output_mesh << "f " << 1 + vertex_count << " " << 2 + vertex_count << " " << 3 + vertex_count << "\n";
			vertex_count += 3;
		}
		if (point_size == 4)
		{
			output_mesh << "v " << points[0*3] << " " << points[0*3+1] << " " << points[0*3+2] << "\n";
			output_mesh << "v " << points[1*3] << " " << points[1*3+1] << " " << points[1*3+2] << "\n";
			output_mesh << "v " << points[2*3] << " " << points[2*3+1] << " " << points[2*3+2] << "\n";
			output_mesh << "v " << points[3*3] << " " << points[3*3+1] << " " << points[3*3+2] << "\n";
			output_mesh << "f " << 1 + vertex_count << " " << 2 + vertex_count << " " << 3 + vertex_count << "\n";
			output_mesh << "f " << 2 + vertex_count << " " << 3 + vertex_count << " " << 4 + vertex_count << "\n";
			output_mesh << "f " << 1 + vertex_count << " " << 2 + vertex_count << " " << 4 + vertex_count << "\n";
			output_mesh << "f " << 3 + vertex_count << " " << 3 + vertex_count << " " << 4 + vertex_count << "\n";
			vertex_count += 4;
		}
	}

	void MCMT::save_triangle_soup(std::string file_name){
		std::ofstream output_mesh(file_name);
		int vertex_count = 0;
		for(index_t t =0; t < delaunay_->nb_finite_cells(); t++)
		{
			std::vector<int> tri2v;
			for (index_t lv = 0; lv < 4; ++lv)
			{
				int v = delaunay_->cell_vertex(t, lv);
				if (v != -1)
					tri2v.push_back(int(v));
			}
			std::vector<double> intersection_points;
			if (point_values_[tri2v[0]] * point_values_[tri2v[1]] < 0)
			{
				std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[0] * 3, point_positions_.data() + tri2v[1] * 3, point_values_[tri2v[0]], point_values_[tri2v[1]]);
				intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
			}
			if (point_values_[tri2v[0]] * point_values_[tri2v[2]] < 0)
			{   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[0] * 3, point_positions_.data() + tri2v[2] * 3, point_values_[tri2v[0]], point_values_[tri2v[2]]);
				intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
			}
			if (point_values_[tri2v[0]] * point_values_[tri2v[3]] < 0)
			{   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[0] * 3, point_positions_.data() + tri2v[3] * 3, point_values_[tri2v[0]], point_values_[tri2v[3]]);
				intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
			}
			if (point_values_[tri2v[1]] * point_values_[tri2v[2]] < 0)
			{   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[1] * 3, point_positions_.data() + tri2v[2] * 3, point_values_[tri2v[1]], point_values_[tri2v[2]]);
				intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
			}
			if (point_values_[tri2v[1]] * point_values_[tri2v[3]] < 0)
			{   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[1] * 3, point_positions_.data() + tri2v[3] * 3, point_values_[tri2v[1]], point_values_[tri2v[3]]);
				intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
			}
			if (point_values_[tri2v[2]] * point_values_[tri2v[3]] < 0)
			{   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[2] * 3, point_positions_.data() + tri2v[3] * 3, point_values_[tri2v[2]], point_values_[tri2v[3]]);
				intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
			}
			if (intersection_points.size() != 0)
			{ 
				save_face(output_mesh, intersection_points, vertex_count);
			}
		} 

	}

	std::vector<double> MCMT::get_mid_points()
	{
		tbb::concurrent_vector<std::vector<double>> sample_points;
		tbb::parallel_for(tbb::blocked_range<int>(0, delaunay_->nb_finite_cells()),
						  [&](tbb::blocked_range<int> ti)
						  {
						  for(index_t t = ti.begin(); t < ti.end(); t++)
						  {
							  std::vector<int> tri2v;
							  for (index_t lv = 0; lv < 4; ++lv)
							  {
								  int v = delaunay_->cell_vertex(t, lv);
								  if (v != -1)
									  tri2v.push_back(int(v));
							  }
							  std::vector<double> intersection_points;
							  if (point_values_[tri2v[0]] * point_values_[tri2v[1]] < 0)
							  {
                                  std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[0] * 3, point_positions_.data() + tri2v[1] * 3, point_values_[tri2v[0]], point_values_[tri2v[1]]);
								  intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
							  }
							  if (point_values_[tri2v[0]] * point_values_[tri2v[2]] < 0)
							  {   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[0] * 3, point_positions_.data() + tri2v[2] * 3, point_values_[tri2v[0]], point_values_[tri2v[2]]);
								  intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
							  }
							  if (point_values_[tri2v[0]] * point_values_[tri2v[3]] < 0)
							  {   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[0] * 3, point_positions_.data() + tri2v[3] * 3, point_values_[tri2v[0]], point_values_[tri2v[3]]);
								  intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
							  }
							  if (point_values_[tri2v[1]] * point_values_[tri2v[2]] < 0)
							  {   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[1] * 3, point_positions_.data() + tri2v[2] * 3, point_values_[tri2v[1]], point_values_[tri2v[2]]);
								  intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
							  }
							  if (point_values_[tri2v[1]] * point_values_[tri2v[3]] < 0)
							  {   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[1] * 3, point_positions_.data() + tri2v[3] * 3, point_values_[tri2v[1]], point_values_[tri2v[3]]);
								  intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
							  }
							  if (point_values_[tri2v[2]] * point_values_[tri2v[3]] < 0)
							  {   std::vector<double> interpolated_point = interpolate(point_positions_.data() + tri2v[2] * 3, point_positions_.data() + tri2v[3] * 3, point_values_[tri2v[2]], point_values_[tri2v[3]]);
								  intersection_points.insert(intersection_points.end(), interpolated_point.begin(), interpolated_point.end());
							  }
							  if (intersection_points.size() != 0)
							  { 
								std::vector<double> mid_point = compute_face_mid_point(intersection_points.size()/3, intersection_points);
							  	sample_points.push_back(mid_point);
							  }
						  } });

		std::vector<double> new_points;
		new_points.reserve(sample_points.size() * 3);
		for (int i = 0; i < sample_points.size(); i++)
		{
			new_points.insert(new_points.end(), sample_points[i].begin(), sample_points[i].end());
		}
		return new_points;
	}

	void MCMT::get_cell(index_t v, ConvexCell &C)
	{
		PeriodicDelaunay3d::IncidentTetrahedra W_;
		delaunay_->copy_Laguerre_cell_from_Delaunay(v, C, W_);
		if (!periodic_)
		{
			C.clip_by_plane(vec4(1.0, 0.0, 0.0, 0.0));
			C.clip_by_plane(vec4(-1.0, 0.0, 0.0, 1.0));
			C.clip_by_plane(vec4(0.0, 1.0, 0.0, 0.0));
			C.clip_by_plane(vec4(0.0, -1.0, 0.0, 1.0));
			C.clip_by_plane(vec4(0.0, 0.0, 1.0, 0.0));
			C.clip_by_plane(vec4(0.0, 0.0, -1.0, 1.0));
		}
		C.compute_geometry();
	}

	std::vector<double> MCMT::lloyd_relaxation(double *relaxed_point_positions, int num_points, int num_iter)
	{
		int current_num_points = point_positions_.size() / 3;
		// copy point_positions_ to new_points
		std::vector<double> new_points;
		new_points.reserve(point_positions_.size() + num_points * 3);
		for (int i = 0; i < point_positions_.size(); i++)
		{
			new_points.push_back(point_positions_[i]);
		}
		for (int i = 0; i < num_points * 3; i++)
		{
			new_points.push_back(relaxed_point_positions[i]);
		}

		ConvexCell C;

		for (int i = 0; i < num_iter; i++)
		{
			delete delaunay_;
			delaunay_ = new PeriodicDelaunay3d(periodic_, 1.0);
			if (!periodic_)
			{
				delaunay_->set_keeps_infinite(true);
			}
			delaunay_->set_vertices(new_points.size() / 3, new_points.data());
			delaunay_->compute();

			std::cout << "Size of points: " << delaunay_->nb_finite_cells() << std::endl;
			for (index_t v = current_num_points; v < new_points.size() / 3; v++)
			{
				get_cell(v, C);
				vec3 g = C.barycenter();
				new_points[3 * v] = g.x;
				new_points[3 * v + 1] = g.y;
				new_points[3 * v + 2] = g.z;
			}
			for (index_t v = current_num_points; v < new_points.size(); v++)
			{
				if (new_points[v] < 0.0)
				{
					new_points[v] += 1.0;
				}
				if (new_points[v] > 1.0)
				{
					new_points[v] -= 1.0;
				}
			}
		}
		std::vector<double> points_vec;
		points_vec.reserve(num_points * 3);
		for (int i = current_num_points*3; i < new_points.size(); i++)
		{
			points_vec.push_back(new_points[i]);
		}
		return points_vec;
	}

	void MCMT::output_grid_points(std::string filename){
		std::ofstream outfile(filename);
		for(int i=0; i<point_positions_.size()/3; i++){
			outfile << "v " << point_positions_[i*3] << " " << point_positions_[i*3+1] << " " << point_positions_[i*3+2] << std::endl;
		}

	}

}
