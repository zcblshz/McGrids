
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <math.h>
#include <random>
#include <algorithm>

// CGAL headers
#define CGAL_LINKED_WITH_TBB
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <tbb/concurrent_unordered_map.h>
struct VertexInfo {
	int point_index;
	double point_value;
	double point_density;
	VertexInfo() : point_index(0), point_value(0.0), point_density(0.0) {}
	VertexInfo(int index, double value, double density)
		: point_index(index), point_value(value), point_density(density) {}
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
typedef CGAL::Triangulation_data_structure_3<
	Vb,
	CGAL::Delaunay_triangulation_cell_base_3<K>,
	CGAL::Parallel_tag>
	Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Finite_cells_iterator Finite_cells_iterator;
typedef Delaunay::Finite_vertices_iterator Finite_vertices_iterator;
typedef Delaunay::Facet Facet;
namespace GEO
{

	class MCMT
	{
	public:
		MCMT();
		~MCMT();
		int add_points(const std::vector<Point> &points, const std::vector<double> &point_values);
		void export_grid_off(const std::string &filename);
		void export_surface_obj(const std::string &filename);

		// void clear();

		// void add_points(int num_points, double *point_positions, double *point_values);
		// void add_mid_points(int num_points, double *point_positions, double *point_values);
		// std::vector<double> get_mid_points();
		// std::vector<double> get_grid_points();
		// std::vector<int> get_grids();
		// std::vector<double> sample_points_rejection(int num_samples, double min_value, double max_value);
		// // std::vector<double> sample_points(int num_samples);
		// std::vector<double> lloyd_relaxation(double *point_positions, int num_points, int num_iter);
		// void output_grid_points(std::string filename);
		// void save_triangle_mesh(std::string filename);
		// void save_grid_mesh(std::string filename, float x_clip_plane);

		// std::vector<double> sample_points_voronoi(const int num_points);

	private:
		Delaunay *delaunay_;
		// PeriodicDelaunay3d::IncidentTetrahedra W_;
		// bool periodic_ = false;
		// double max_bound = 0;
		// double min_bound = 0;
		// int num_point_visited_ = 0;
		// std::vector<double> point_positions_;
		// std::vector<double> point_values_;
		// std::vector<double> point_errors_;
		// std::vector<double> point_volumes_;
		// std::vector<bool> volume_changed_;

		// std::vector<double> sample_tet(std::vector<double> point_positions);
		// std::vector<double> compute_tet_error();
		// std::vector<double> sample_polytope(int vertex_index);
		// std::vector<double> compute_voronoi_error();

		// double tetrahedronVolume(const std::vector<double> &coordinates);
		// void save_face(std::ofstream &output_mesh, const std::vector<double> &points, int &vertex_count);
		// void get_cell(index_t v, ConvexCell &C, PeriodicDelaunay3d::IncidentTetrahedra& W);

		// std::vector<double> compute_face_mid_point(int num_points, const std::vector<double> &points);
		Point interpolate(Vertex_handle p1, Vertex_handle p2);
		std::vector<double> compute_point_bbox(const std::vector<Point>& points);
		double compute_point_density(double point_value);
	};
}