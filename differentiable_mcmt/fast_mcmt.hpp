
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
	bool visited = false;
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
		void clear();
		int add_points(const std::vector<Point> &points, const std::vector<double> &point_values);
		std::vector<Point> get_mid_points();
		void export_grid_off(const std::string &filename);
		void export_surface_obj(const std::string &filename);
	private:
		Delaunay *delaunay_;
		int num_point_visited_ = 0;
		std::vector<std::pair<int, int>> configurations_;
		Point interpolate(Vertex_handle p1, Vertex_handle p2);
		std::vector<double> compute_point_bbox(const std::vector<Point>& points);
		double compute_point_density(double point_value);
	};
}