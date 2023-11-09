#include "fast_cvt3d.h"

/*
float ScTP_(Vector_3 a, Vector_3 b, Vector_3 c)
{
	// computes scalar triple product
	return a * CGAL::cross_product(b, c);
}

std::vector<double> computeBC_(Point v0, Point v1, Point v2, Point v3, Point p) {

	Vector_3 vap = Point_3(p.x(), p.y(), p.z()) - Point_3(v0.x(), v0.y(), v0.z());
	Vector_3 vbp = Point_3(p.x(), p.y(), p.z()) - Point_3(v1.x(), v1.y(), v1.z());

	Vector_3 vab = Point_3(v1.x(), v1.y(), v1.z()) - Point_3(v0.x(), v0.y(), v0.z());
	Vector_3 vac = Point_3(v2.x(), v2.y(), v2.z()) - Point_3(v0.x(), v0.y(), v0.z());
	Vector_3 vad = Point_3(v3.x(), v3.y(), v3.z()) - Point_3(v0.x(), v0.y(), v0.z());

	Vector_3 vbc = Point_3(v2.x(), v2.y(), v2.z()) - Point_3(v1.x(), v1.y(), v1.z());
	Vector_3 vbd = Point_3(v3.x(), v3.y(), v3.z()) - Point_3(v1.x(), v1.y(), v1.z());

	// ScTP computes the scalar triple product
	float va6 = ScTP_(vbp, vbd, vbc);
	float vb6 = ScTP_(vap, vac, vad);
	float vc6 = ScTP_(vap, vad, vab);
	float vd6 = ScTP_(vap, vab, vac);
	float v6 = 1 / ScTP_(vab, vac, vad);

	return std::vector<double>{ va6*v6, vb6*v6, vc6*v6, vd6*v6 };
}

Point calculateCentroid(Point A, Point B, Point C, Point D) {

	float x = (A.x() + B.x() + C.x() + D.x()) / 4;
	float y = (A.y() + B.y() + C.y() + D.y()) / 4;
	float z = (A.z() + B.z() + C.z() + D.z()) / 4;
	return Point(x, y, z);
}

double calculateTetrahedronVolume(Point A, Point B, Point C, Point D) {
	Vector_3 diff1 = { A.x() - D.x(), A.y() - D.y(), A.z() - D.z() };
	Vector_3 diff2 = { B.x() - D.x(), B.y() - D.y(), B.z() - D.z() };
	Vector_3 diff3 = { C.x() - D.x(), C.y() - D.y(), C.z() - D.z() };

	double volume = std::abs(diff1 * CGAL::cross_product(diff2, diff3)) / 6;
	return volume;
}
*/

namespace GEO
{
	Fast_CVT::Fast_CVT()
	{
		// Constructor implementation
		// delaunay_ = Delaunay::create(dimension_, "PDEL");
		delaunay_ = new PeriodicDelaunay3d(periodic_, 1.0);
		if (!periodic_)
		{
			delaunay_->set_keeps_infinite(true);
		}
	}

	Fast_CVT::~Fast_CVT()
	{
		// Destructor implementation
	}

	void Fast_CVT::set_points(index_t nb_total_points, index_t nb_constrain_points, const double *points)
	{
		// clear points_
		// points_.resize(dimension_ * nb_total_points);
		std::cout << "=============" << std::endl;
		std::cout << "points_ size: " << points_.size()/3 << std::endl;
		for (index_t i = 0; i < nb_total_points; i++)
		{
			points_.push_back(points[i*3]);
			points_.push_back(points[i*3+1]);
			points_.push_back(points[i*3+2]);

			if (i < nb_constrain_points)
				lock_point(i);
		}
		delaunay_->set_vertices(points_.size()/3, points_.data());
		delaunay_->compute();
		std::cout << "delaunay_ finished" << std::endl;

	}

	void Fast_CVT::set_errors(int nb_errors, const float *errors)
	{
		for (int i = 0; i < nb_errors; i++)
		{
			point_error_.push_back((double)errors[i]);
		}
		delaunay_->set_weights(point_error_.data());
	}

	void Fast_CVT::set_sdfs(int nb_errors, const float *sdfs)
	{
		std::cout << "sdf size: " << point_sdfs_.size() << std::endl;
		
		for (int i = 0; i < nb_errors; i++)
		{
			point_sdfs_.push_back((double)sdfs[i]);
			if (sdfs[i] < 0)
				point_signs_.push_back(-1);
			else
				point_signs_.push_back(1);
		}
	}

	void Fast_CVT::get_cell(index_t v, ConvexCell &C)
	{
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

	void Fast_CVT::Lloyd_iterations(index_t nb_iter)
	{

		index_t nb_points = index_t(points_.size() / dimension_);

		vector<double> new_points(points_.size());
		ConvexCell C;

		cur_iter_ = 0;
		nb_iter_ = nb_iter;

		for (index_t cur_iter_ = 0; cur_iter_ < nb_iter; cur_iter_++)
		{
			for (index_t v = 0; v < points_.size() / 3; ++v)
			{
				if (!point_is_locked(v))
				{
					get_cell(v, C);
					vec3 g = C.barycenter();
					new_points[3 * v] = g.x;
					new_points[3 * v + 1] = g.y;
					new_points[3 * v + 2] = g.z;
				}
				else
				{
					new_points[3 * v] = points_[3 * v];
					new_points[3 * v + 1] = points_[3 * v + 1];
					new_points[3 * v + 2] = points_[3 * v + 2];
				}
			}
			// In periodic mode, points may escape out of
			// the domain. Relocate them in [0,1]^3
			for (index_t i = 0; i < new_points.size(); ++i)
			{
				if (new_points[i] < 0.0)
				{
					new_points[i] += 1.0;
				}
				if (new_points[i] > 1.0)
				{
					new_points[i] -= 1.0;
				}
			}
			points_.swap(new_points);
			delaunay_->set_vertices(points_.size() / 3, points_.data());
			delaunay_->compute();
			cur_iter_++;
		}
	}

	void Fast_CVT::output_delaunay(std::string filename)
	{
		tetra_vertices.clear();
		tetra_indices.clear();

		std::ofstream outfile(filename); // create a file named "example.txt"

		if (outfile.is_open())
		{

			const double *point = delaunay_->vertices_ptr();
			for (index_t t = 0; t < delaunay_->nb_vertices(); ++t)
			{
				double x = point[t * 3];
				double y = point[t * 3 + 1];
				double z = point[t * 3 + 2];
				tetra_vertices.push_back(std::vector<double>{x, y, z});

				outfile << "v " << x << " " << y << " " << z << " \n";
			}

			for (index_t t = 0; t < delaunay_->nb_cells(); t++)
			{
				std::vector<int> tri2v;
				for (index_t lv = 0; lv < 4; ++lv)
				{
					int v = delaunay_->cell_vertex(t, lv);
					if (v != -1)
					{
						tri2v.push_back(int(v));
					}
				}
				tetra_indices.push_back(tri2v);
				if (tri2v.size() != 4)
				{
					continue;
				}
				outfile << "f " << tri2v[0] + 1 << " " << tri2v[2] + 1 << " " << tri2v[1] + 1 << " \n";
				outfile << "f " << tri2v[0] + 1 << " " << tri2v[3] + 1 << " " << tri2v[2] + 1 << " \n";
				outfile << "f " << tri2v[0] + 1 << " " << tri2v[1] + 1 << " " << tri2v[3] + 1 << " \n";
				outfile << "f " << tri2v[1] + 1 << " " << tri2v[2] + 1 << " " << tri2v[3] + 1 << " \n";
			}
		}
	}

	std::vector<vec3> Fast_CVT::near_surface_sampling(int nb_points)
	{
		// calculate tetrahedron errors
		std::vector<std::vector<int>> tet_lists;
		std::vector<float> tet_errors;
		float total_error = 0.0;
		for (index_t t = 0; t < delaunay_->nb_cells(); t++)
		{
			std::vector<int> tri2v;
			float error = 0.0;
			for (index_t lv = 0; lv < 4; ++lv)
			{
				int v = delaunay_->cell_vertex(t, lv);
				if (v != -1)
				{
					tri2v.push_back(int(v));
					error += point_error_[int(v)];
				}
			}
			if (tri2v.size() != 4)
			{
				continue;
			}
			tet_lists.push_back(tri2v);
			tet_errors.push_back(error / 4.);
			total_error += error / 4.;
		}

		// normalize tetrahedron errors
		float curr_error = 0.0;
		for (int i = 0; i < tet_lists.size(); i++)
		{
			tet_errors[i] = curr_error + tet_errors[i] / total_error;
			curr_error = tet_errors[i];
		}

		// randomly sampling
		std::vector<vec3> samples;
		for (int i = 0; i < nb_points; i++)
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<double> dis(0.0, 1.0);

			double rand = dis(gen);

			std::vector<float>::iterator upper;
			upper = std::upper_bound(tet_errors.begin(), tet_errors.end(), rand);
			int idx = upper - tet_errors.begin();

			std::vector<int> tri2v = tet_lists[idx];
			if (tri2v.size() != 4)
				continue;
			vec3 sample = sample_tet(tri2v);
			samples.push_back(sample);
		}

		return samples;
	}

	vec3 Fast_CVT::sample_tet(std::vector<int> tri2v)
	{
		const vec3 A = delaunay_->vertex(tri2v[0]);
		const vec3 B = delaunay_->vertex(tri2v[1]);
		const vec3 C = delaunay_->vertex(tri2v[2]);
		const vec3 D = delaunay_->vertex(tri2v[3]);

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dis(0.0, 1.0);

		double s = dis(gen);
		double t = dis(gen);
		double u = dis(gen);

		if (s + t > 1)
		{
			s = 1 - s;
			t = 1 - t;
		}

		if (t + u > 1)
		{
			double tmp = u;
			u = 1 - s - t;
			t = 1 - tmp;
		}
		else if (s + t + u > 1)
		{
			double tmp = u;
			u = s + t + u - 1;
			s = 1 - t - tmp;
		}

		double w = 1 - s - t - u;
		vec3 sample(s * A.x + t * B.x + u * C.x + w * D.x,
					s * A.y + t * B.y + u * C.y + w * D.y,
					s * A.z + t * B.z + u * C.z + w * D.z);

		return sample;
	}

	// std::vector<int> Fast_CVT::get_cell_indices(int nb_points)
	// {
	// 	float total_error = 0.0;
	// 	for (int i = 0; i < points_.size(); i++)
	// 	{
	// 		total_error += point_error_[i];
	// 	}
	// 	vector<float> normalized_point_error_;
	// 	normalized_point_error_.resize(points_.size());
	// 	for (int i = 0; i < points_.size(); i++)
	// 	{
	// 		normalized_point_error_[i] = point_error_[i] / total_error;
	// 	}

	// 	std::vector<int> sample_indices;
	// 	for (int i = 0; i < nb_points; i++)
	// 	{
	// 		float rand = Numeric::random_float64();

	// 		std::vector<float>::iterator upper;
	// 		upper = std::upper_bound(normalized_point_error_.begin(), normalized_point_error_.end(), rand);
	// 		int idx = upper - normalized_point_error_.begin();
	// 		sample_indices.push_back(idx);
	// 	}

	// 	return sample_indices;
	// }

	std::vector<vec3> Fast_CVT::sample_cells(int *cell_indices, int nb_cell)
	{
		std::vector<vec3> sample_points;
		ConvexCell C;
		std::ofstream outfile("cell.obj"); // create a file named "example.txt"

		for (int i = 0; i < nb_cell; i++)
		{
			std::vector<vec3> vertices;
			get_cell(cell_indices[i], C);

			for (index_t j = 0; j < C.nb_v(); j++)
			{

				vec3 v = C.stored_triangle_point(j);

				if ((v.x == 0) && (v.y == 0) && (v.z == 0))
					continue;

				vertices.push_back(v);
				outfile << "v " << v.x << " " << v.y << " " << v.z << "\n";
			}
			std::vector<float> weights;
			float sum = 0;
			for (int w = 0; w < vertices.size(); w++)
			{
				float weight = Numeric::random_float64();
				sum += weight;
				weights.push_back(weight);
			}
			vec3 sample_point(vertices[0].x * weights[0] / sum, vertices[0].y * weights[0] / sum, vertices[0].z * weights[0] / sum);
			for (int j = 1; j < weights.size(); j++)
			{
				sample_point.x += vertices[j].x * weights[j] / sum;
				sample_point.y += vertices[j].y * weights[j] / sum;
				sample_point.z += vertices[j].z * weights[j] / sum;
			}
			sample_points.push_back(sample_point);
		}
		return sample_points;
	}
	vec3 Fast_CVT::interpolate(double* ip1, double* ip2, double sd1, double sd2)
	{
		vec3 p1(ip1[0], ip1[1], ip1[2]);
		vec3 p2(ip2[0], ip2[1], ip2[2]);
		double t = sd1 / ((sd1 - sd2) + 1e-6);
		return vec3(p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y), p1.z + t * (p2.z - p1.z));
	}

	vec3 Fast_CVT::compute_mid_point(const std::vector<vec3> &points)
	{
		// std::cout << "points: " << points.size() << std::endl;
		// std::cout << "point_0" << points[0].x << " " << points[0].y << " " << points[0].z << std::endl;
		vec3 mid_point(0, 0, 0);
		for (int i=0; i < points.size(); i++)
		{
			mid_point.x += points[i].x;
			mid_point.y += points[i].y;
			mid_point.z += points[i].z;

		}

		mid_point.x /= static_cast<double>(points.size());
		mid_point.y /= static_cast<double>(points.size());
		mid_point.z /= static_cast<double>(points.size());
		return mid_point;
	}

	void Fast_CVT::add_points(const std::vector<vec3>& points, const std::vector<double>&sdfs){
		for (int i=0; i<points.size(); i++){
			points_.push_back(points[i].x);
			points_.push_back(points[i].y);
			points_.push_back(points[i].z);
			point_sdfs_.push_back(sdfs[i]);
			if (sdfs[i] < 0)
				point_signs_.push_back(-1);
			else
				point_signs_.push_back(1);
		}
		delaunay_->set_vertices(points_.size() / 3, points_.data());
		delaunay_->compute();
	

	}

	std::vector<vec3> Fast_CVT::compute_cell_intersection_mid_point()
	{
		std::vector<vec3> new_points;
		for (index_t t = 0; t < delaunay_->nb_finite_cells(); t++)
		{
			std::vector<int> tri2v;
			for (index_t lv = 0; lv < 4; ++lv)
			{
				int v = delaunay_->cell_vertex(t, lv);
				if (v != -1)
					tri2v.push_back(int(v));
			}
			double error = 0;
			std::vector<vec3> intersection_points;
			if (point_signs_[tri2v[0]] * point_signs_[tri2v[1]] < 0)
			{
				intersection_points.push_back(interpolate(points_.data() + tri2v[0]*3, points_.data()+tri2v[1]*3, point_sdfs_[tri2v[0]], point_sdfs_[tri2v[1]]));
			}
			if (point_signs_[tri2v[0]] * point_signs_[tri2v[2]] < 0)
			{
				intersection_points.push_back(interpolate(points_.data() + tri2v[0]*3, points_.data()+tri2v[2]*3, point_sdfs_[tri2v[0]], point_sdfs_[tri2v[2]]));
			}
			if (point_signs_[tri2v[0]] * point_signs_[tri2v[3]] < 0)
			{
				intersection_points.push_back(interpolate(points_.data() + tri2v[0]*3, points_.data()+tri2v[3]*3, point_sdfs_[tri2v[0]], point_sdfs_[tri2v[3]]));
			}
			if (point_signs_[tri2v[1]] * point_signs_[tri2v[2]] < 0)
			{
				intersection_points.push_back(interpolate(points_.data() + tri2v[1]*3, points_.data()+tri2v[2]*3, point_sdfs_[tri2v[1]], point_sdfs_[tri2v[2]]));
			}
			if (point_signs_[tri2v[1]] * point_signs_[tri2v[3]] < 0)
			{
				intersection_points.push_back(interpolate(points_.data() + tri2v[1]*3, points_.data()+tri2v[3]*3, point_sdfs_[tri2v[1]], point_sdfs_[tri2v[3]]));
			}
			if (point_signs_[tri2v[2]] * point_signs_[tri2v[3]] < 0)
			{
				intersection_points.push_back(interpolate(points_.data() + tri2v[2]*3, points_.data()+tri2v[3]*3, point_sdfs_[tri2v[2]], point_sdfs_[tri2v[3]]));
			}
			if (intersection_points.size() == 0)
			{
				continue;
			}
			vec3 mid_point = compute_mid_point(intersection_points);
			new_points.push_back(mid_point);
		}
		return new_points;
	}



	

}
/*
void Fast_CVT::outputClippedTesellation(Delaunay T_delaunay, std::string filename)
{
	//clipped x < 0 tetrahedron

	tetra_vertices.clear();
	tetra_indices.clear();

	std::map<Vertex_handle, int> vh_map;
	std::map<Cell_handle, int> ch_map;

	int vcounter = 0;
	Delaunay::Finite_vertices_iterator vits = T_delaunay.finite_vertices_begin();
	Delaunay::Finite_vertices_iterator vite = T_delaunay.finite_vertices_end();
	for (; vits != vite; ++vits)
	{
		Vertex_handle vh = Vertex_handle(vits);
		Point p = vh->point();

		if (p.x() < 0)
			continue;

		vh_map.insert(std::pair<Vertex_handle, int>(vh, vcounter));
		vcounter++;
		tetra_vertices.push_back(std::vector<float>{float(p.x()), float(p.y()), float(p.z())});
	}

	int ccounter = 0;
	Delaunay::Finite_cells_iterator cits = T_delaunay.finite_cells_begin();
	Delaunay::Finite_cells_iterator cite = T_delaunay.finite_cells_end();
	for (; cits != cite; ++cits)
	{
		Cell cell = Cell(*cits);
		ch_map.insert(std::pair<Cell_handle, int>(Cell_handle(cits), ccounter));
		ccounter++;
		std::vector<int> vindices;
		bool v_flag = true;
		for (int i = 0; i < 4; i++)
		{
			Vertex_handle vh = cell.vertex(i);
			if (vh_map.find(vh) == vh_map.end()) {
				// not found
				v_flag = false;
				break;
			}
			else {
				int vid = vh_map[vh];
				vindices.push_back(vid);
			}

		}
		if(v_flag)
			tetra_indices.push_back(vindices);
	}

	std::ofstream outfile(filename); // create a file named "example.txt"

	if (outfile.is_open()) {
		for (size_t k = 0; k < tetra_vertices.size(); ++k) {
			std::vector<float> site = tetra_vertices[k];
			size_t n = site.size();
			outfile << "v ";
			for (size_t i = 0; i < site.size(); i++) {
				outfile << site[i] << " ";

			}
			outfile << " \n";
		}

		for (size_t j = 0; j < tetra_indices.size(); j++) {
			std::vector<int> vertices_idx = tetra_indices[j];
			size_t n = vertices_idx.size();

			if (n == 4) {
				outfile << "f " << vertices_idx[0] + 1 << " " << vertices_idx[2] + 1 << " " << vertices_idx[1] + 1 << " \n";
				outfile << "f " << vertices_idx[0] + 1 << " " << vertices_idx[3] + 1 << " " << vertices_idx[2] + 1 << " \n";
				outfile << "f " << vertices_idx[0] + 1 << " " << vertices_idx[1] + 1 << " " << vertices_idx[3] + 1 << " \n";
				outfile << "f " << vertices_idx[1] + 1 << " " << vertices_idx[2] + 1 << " " << vertices_idx[3] + 1 << " \n";
			}
		}
	}

}

*/