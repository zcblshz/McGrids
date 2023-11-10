#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "fast_mcmt.hpp"

namespace py = pybind11;

namespace mcmt
{
  GEO::MCMT mcmt = GEO::MCMT();

  void add_points(py::array_t<double> &point_positions, py::array_t<double> &point_values)
  {
    py::buffer_info point_positions_buffer = point_positions.request();
    double *point_positions_prt = (double *)point_positions_buffer.ptr;
    int num_points = point_positions_buffer.shape[0];
    py::buffer_info point_values_buffer = point_values.request();
    double *point_values_prt = (double *)point_values_buffer.ptr;
    int num_values = point_values_buffer.shape[0];
    if (num_points != num_values)
    {
      throw std::runtime_error("Input points and values must have the same number of rows");
    }
    mcmt.add_points(num_points, point_positions_prt, point_values_prt);
  }

  static py::array_t<double> sample_points(int num_points)
  {
    std::vector<double> new_samples = mcmt.sample_points(num_points);
    return py::array_t<double>(new_samples.size(), new_samples.data());
  }

  static py::array_t<double> get_grid_points()
  {
    std::vector<double> grid_points = mcmt.get_grid_points();
    return py::array_t<double>(grid_points.size(), grid_points.data());
  }

  static py::array_t<double> lloyd_relaxation(py::array_t<double> &point_positions, int num_iter)
  {
    py::buffer_info point_positions_buffer = point_positions.request();
    double *point_positions_prt = (double *)point_positions_buffer.ptr;
    int num_points = point_positions_buffer.shape[0];
    std::vector<double> new_samples = mcmt.lloyd_relaxation(point_positions_prt, num_points, num_iter);
    return py::array_t<double>(new_samples.size(), new_samples.data());
  }

  static py::array_t<double> get_mid_points()
  {
    std::vector<double> mid_points = mcmt.get_mid_points();
    return py::array_t<double>(mid_points.size(), mid_points.data());
  }

  void output_triangle_soup(std::string filename)
  {
    mcmt.save_triangle_soup(filename);
  }

  void output_triangle_mesh(std::string filename)
  {
    mcmt.save_triangle_mesh(filename);
  }


  PYBIND11_MODULE(mcmt, m)
  {

    m.def("add_points", &add_points,
          "add points to mcmt");
    m.def("sample_points", &sample_points,
          "monte carlo sampling from mcm");
    m.def("lloyd_relaxation", &lloyd_relaxation,
          "lloyd relaxation");
    m.def("get_grid_points", &get_grid_points,
          "get all grid points");
    m.def("get_mid_points", &get_mid_points,
          "get all mid points");
    m.def("output_triangle_soup", &output_triangle_soup,
          "output triangle soup");
    m.def("output_triangle_mesh", &output_triangle_mesh,
          "output triangle mesh");
  }

}
