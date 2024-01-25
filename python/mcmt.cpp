#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "fast_mcmt.hpp"

namespace py = pybind11;

namespace mcmt
{
  MCMT::MCMT mcmt = MCMT::MCMT();

  void add_points(py::array_t<double> &point_positions, py::array_t<double> &point_values, py::array_t<double> &point_curvatures)
  {
    py::buffer_info point_positions_buffer = point_positions.request();
    double *point_positions_prt = (double *)point_positions_buffer.ptr;
    int num_points = point_positions_buffer.shape[0];
    py::buffer_info point_values_buffer = point_values.request();
    double *point_values_prt = (double *)point_values_buffer.ptr;

    py::buffer_info point_curvatures_buffer = point_curvatures.request();
    double *point_curvatures_prt = (double *)point_curvatures_buffer.ptr;

    int num_values = point_values_buffer.shape[0];
    if (num_points != num_values)
    {
      throw std::runtime_error("Input points and values must have the same number of rows");
    }
    std::vector<Point> points;
    std::vector<double> point_values_vec;
    std::vector<double> point_curvatures_vec;

    for (int i = 0; i < num_points; i++)
    {
      points.push_back(Point(point_positions_prt[3 * i], point_positions_prt[3 * i + 1], point_positions_prt[3 * i + 2]));
      point_values_vec.push_back(point_values_prt[i]);
      point_curvatures_vec.push_back(point_curvatures_prt[i]);
    }
    mcmt.add_points(points, point_values_vec, point_curvatures_vec);
  }
  
  static py::array_t<double> sample_points_tetrahedron(int num_points)
  {
    std::vector<Point> new_samples = mcmt.sample_tetrahedron(num_points);
    std::vector<double> new_samples_vec;
    for (int i = 0; i < new_samples.size(); i++)
    {
      new_samples_vec.push_back(new_samples[i].x());
      new_samples_vec.push_back(new_samples[i].y());
      new_samples_vec.push_back(new_samples[i].z());
    }
    return py::array_t<double>(new_samples_vec.size(), new_samples_vec.data());
  }

  static py::array_t<double> get_mid_points()
  {
    std::vector<Point> mid_points = mcmt.get_mid_points();
    std::vector<double> mid_points_vec;
    for (int i = 0; i < mid_points.size(); i++)
    {
      mid_points_vec.push_back(mid_points[i].x());
      mid_points_vec.push_back(mid_points[i].y());
      mid_points_vec.push_back(mid_points[i].z());
    }
    return py::array_t<double>(mid_points_vec.size(), mid_points_vec.data());
  }

  void output_triangle_mesh(std::string filename)
  {
    mcmt.export_surface_obj(filename);
  }

  void output_grid_mesh(std::string filename)
  {
    mcmt.export_grid_off(filename);
  }

  void clear_mcmt()
  {
    mcmt.clear();
  }

  PYBIND11_MODULE(mcmt, m)
  {

    m.def("add_points", &add_points,
          "add points to mcmt");
    m.def("sample_points_tetrahedron", &sample_points_tetrahedron,
          "monte carlo sampling from mcm");
    m.def("get_mid_points", &get_mid_points,
          "get all mid points");
    m.def("output_triangle_mesh", &output_triangle_mesh,
          "output triangle mesh");
    m.def("output_grid_mesh", &output_grid_mesh,
          "output grid mesh");
    m.def("clear", &clear_mcmt,
          "clear mcmt");
  }

}
