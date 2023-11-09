
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <math.h>
#include <random>
#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/periodic_delaunay_3d.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/generic_RVD.h>
#include <geogram/voronoi/integration_simplex.h>
#include <geogram/voronoi/convex_cell.h>
#include <algorithm>

namespace GEO {

	class Fast_CVT
	{
	public:

		bool tbb = false;

		Fast_CVT();
		~Fast_CVT();

		void set_points(index_t nb_points, index_t nb_constrain_points, const double* points);
		void set_errors(int nb_errors, const float *errors);

		// cvt relaxation for point sets
		void Lloyd_iterations(index_t nb_iter);

		// First step near surface sampling
		std::vector<vec3> near_surface_sampling(int nb_points);

		// Second step intersection sampling
		std::vector<vec3> intersection_sampling(int nb_points);

		// IO functions
		void output_delaunay(std::string filename);


	private:

		PeriodicDelaunay3d* delaunay_;
		PeriodicDelaunay3d::IncidentTetrahedra W_;

		vector<double> points_;
		std::vector<double> point_error_;
		vector<bool> point_is_locked_;

		index_t cur_iter_;
		index_t nb_iter_;

		bool periodic_ = false;
		coord_index_t dimension_ = 3;

		std::vector<std::vector<double>> tetra_vertices;
		std::vector<std::vector<int>> tetra_indices;

		index_t nb_points() const {
			return index_t(points_.size() / dimension_);
		}

		void get_cell(index_t v, ConvexCell& C);

		// cvt related functions
		bool point_is_locked(index_t i) const {
			return point_is_locked_.size() != 0 && point_is_locked_[i];
		}
		void lock_point(index_t i) {
			if (point_is_locked_.size() != nb_points()) {
				point_is_locked_.resize(nb_points(), false);
			}
			point_is_locked_[i] = true;
		}
		void unlock_point(index_t i) {
			point_is_locked_[i] = false;
		}
		void unlock_all_points() {
			point_is_locked_.clear();
		}


		// First step private functions
		std::vector<int>  get_cell_indices(int nb_points);
		std::vector<vec3> sample_cells(int* cell_indices, int nb_cell);
		std::vector<int> Fast_CVT::get_tets(int nb_points);
		vec3 sample_tet(std::vector<int> tri2v);

	};
}