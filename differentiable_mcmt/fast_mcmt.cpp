#include "fast_mcmt.hpp"
#include <tbb/tbb.h>

namespace GEO
{
    MCMT::MCMT()
    {
        delaunay_ = nullptr;
        configurations_.clear();
        configurations_.push_back(std::make_pair(0, 1));
        configurations_.push_back(std::make_pair(0, 2));
        configurations_.push_back(std::make_pair(0, 3));
        configurations_.push_back(std::make_pair(1, 2));
        configurations_.push_back(std::make_pair(1, 3));
        configurations_.push_back(std::make_pair(2, 3));
    }

    MCMT::~MCMT()
    {
        delete delaunay_;
    }

    int MCMT::add_points(const std::vector<Point> &points, const std::vector<double> &point_values)
    {
        int current_num_vertices = 0;
        if (delaunay_ != nullptr)
            current_num_vertices = delaunay_->number_of_vertices();
        std::vector<double> bbox = compute_point_bbox(points);
        std::vector<std::pair<Point, VertexInfo>> points_with_info;
        for (int i = 0; i < points.size(); i++)
        {
            points_with_info.push_back(std::make_pair(points[i], VertexInfo(i + current_num_vertices, point_values[i], compute_point_density(point_values[i]))));
        }
        // create a lock datastructure for parallel insertion
        Delaunay::Lock_data_structure locking_ds(CGAL::Bbox_3(bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]), 50);
        if (delaunay_ == nullptr)
            delaunay_ = new Delaunay(points_with_info.begin(), points_with_info.end(), &locking_ds);
        else
            delaunay_->insert(points_with_info.begin(), points_with_info.end(), &locking_ds);

        assert(delaunay_->is_valid());
        return delaunay_->number_of_vertices();
    }

    void MCMT::clear()
    {
        delete delaunay_;
    }

    std::vector<Point> MCMT::sample_tetrahedron(int num_points)
    {
        tbb::concurrent_vector<double> tetrahedron_volume;
        tetrahedron_volume.resize(delaunay_->number_of_finite_cells());
        std::vector<Finite_cells_iterator> cells;
        for (auto it = delaunay_->finite_cells_begin(); it != delaunay_->finite_cells_end(); ++it)
        {
            cells.push_back(it);
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, delaunay_->number_of_finite_cells()),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t i = r.begin(); i != r.end(); ++i)
                {
                    const Finite_cells_iterator& cit = cells[i];
                    double density_sum = cit->vertex(0)->info().point_density + cit->vertex(1)->info().point_density + cit->vertex(2)->info().point_density + cit->vertex(3)->info().point_density;
                    Point density = Point(cit->vertex(0)->info().point_density, cit->vertex(1)->info().point_density, cit->vertex(2)->info().point_density, cit->vertex(3)->info().point_density);
                    tetrahedron_volume[i] = delaunay_->tetrahedron(cit).volume() * density_sum;
                }
            });

        std::cout << "012" << std::endl;

        // normalize the volume
        double total_volume = tbb::parallel_reduce(tbb::blocked_range<tbb::concurrent_vector<double>::iterator>(tetrahedron_volume.begin(), tetrahedron_volume.end()), 0.0,
                                        [](tbb::blocked_range<tbb::concurrent_vector<double>::iterator> const& range, double init) {
                                            return std::accumulate(range.begin(), range.end(), init);
                                        }, std::plus<double>()
                                        );

        tbb::parallel_for(tbb::blocked_range<size_t>(0, tetrahedron_volume.size()), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                tetrahedron_volume[i] /= total_volume;
            }
        });

        // sample tetrahedron
        std::random_device rd;
        std::mt19937 gen(rd());
        std::discrete_distribution<int> volume_distribution(tetrahedron_volume.begin(), tetrahedron_volume.end());
        tbb::concurrent_vector<Point> sampled_points;
        
        tbb::parallel_for(0, num_points, [&](int n) {
            int cell_index = volume_distribution(gen);
            Point sample_point = sample_tetrahedron(cells[cell_index]);
            sampled_points.push_back(sample_point);
        });

        // convert concurrent_vector to vector
        std::vector<Point> sampled_points_vec(sampled_points.begin(), sampled_points.end());
        return sampled_points_vec;
    }

    std::vector<Point> MCMT::get_mid_points()
    {
        tbb::concurrent_vector<Cell_handle> new_cells;

        std::vector<Finite_cells_iterator> cells;
        for (auto it = delaunay_->finite_cells_begin(); it != delaunay_->finite_cells_end(); ++it)
        {
            cells.push_back(it);
        }
        tbb::parallel_for_each(
            cells.begin(),
            cells.end(),
            [&](const Finite_cells_iterator &cit)
            {
                bool skip_mid_point = false;
                // if all vertices are visited, skip
                if (cit->vertex(0)->info().visited && cit->vertex(1)->info().visited && cit->vertex(2)->info().visited && cit->vertex(3)->info().visited)
                {
                    skip_mid_point = true;
                    return;
                }
                // if volume is too small, skip
                double volume = delaunay_->tetrahedron(cit).volume();
                if (volume < 1e-9)
                {
                    skip_mid_point = true;
                    return;

                }

                // if all vertices are positive or negative, skip
                unsigned char index = 0;
                for (int lv = 0; lv < 4; ++lv)
                {
                    if (cit->vertex(lv)->info().point_value < 0)
                        index |= (1 << lv);
                }
                if (index == 0x00 || index == 0x0F)
                {
                    skip_mid_point = true;
                    return;
                }
                // if all pass, compute the mid point
                if (!skip_mid_point)
                    new_cells.push_back(cit);
            });

        tbb::concurrent_vector<Point> mid_points;

        tbb::parallel_for_each(
            new_cells.begin(),
            new_cells.end(),
            [&](const Cell_handle &ch)
            {
                // set all vertex as visited
                for (int lv = 0; lv < 4; ++lv)
                {
                    ch->vertex(lv)->info().visited = true;
                }

                double x_sum = 0, y_sum = 0, z_sum = 0;
                int num_intersection = 0;
                for (std::pair<int, int> config : configurations_)
                {
                    int v1 = config.first;
                    int v2 = config.second;
                    if (ch->vertex(v1)->info().point_value * ch->vertex(v2)->info().point_value < 0)
                    {
                        Point intersection_point = interpolate(ch->vertex(v1), ch->vertex(v2));
                        x_sum += intersection_point.x();
                        y_sum += intersection_point.y();
                        z_sum += intersection_point.z();
                        num_intersection++;
                    }
                }
                mid_points.push_back(Point(x_sum / num_intersection, y_sum / num_intersection, z_sum / num_intersection));
            });
        // convert concurrent_vector to vector
        std::vector<Point> mid_points_vec(mid_points.begin(), mid_points.end());
        return mid_points_vec;
    }

    Point MCMT::sample_tetrahedron(Cell_handle cell)
    {
        Point p0 = cell->vertex(0)->point();
        Point p1 = cell->vertex(1)->point();
        Point p2 = cell->vertex(2)->point();
        Point p3 = cell->vertex(3)->point();
        double density0 = cell->vertex(0)->info().point_density;
        double density1 = cell->vertex(1)->info().point_density;
        double density2 = cell->vertex(2)->info().point_density;
        double density3 = cell->vertex(3)->info().point_density;
        double density_sum = density0 + density1 + density2 + density3;
        density0 /= density_sum;
        density1 /= density_sum;
        density2 /= density_sum;
        density3 /= density_sum;
        // generate 3 random number ranging 0-1
        double s = (double)rand() / RAND_MAX;
        double t = (double)rand() / RAND_MAX;
        double u = (double)rand() / RAND_MAX;

        if (s + t > 1.0)
        {
            s = 1.0 - s;
            t = 1.0 - t;
        }
        if (t + u > 1.0)
        {
            double tmp = u;
            u = 1.0 - s - t;
            t = 1.0 - tmp;
        }
        else if (s + t + u > 1.0)
        {
            double tmp = u;
            u = s + t + u - 1.0;
            s = 1 - t - tmp;
        }
        double a = 1 - s - t - u;
        a *= density0;
        s *= density1;
        t *= density2;
        u *= density3;

        double sum = a + s + t + u;
        a /= sum;
        s /= sum;
        t /= sum;
        u /= sum;

        double x = a * p0.x() + s * p1.x() + t * p2.x() + u * p3.x();
        double y = a * p0.y() + s * p1.y() + t * p2.y() + u * p3.y();
        double z = a * p0.z() + s * p1.z() + t * p2.z() + u * p3.z();

        return Point(x, y, z);
    }

    void MCMT::export_grid_off(const std::string &filename)
    {
        std::ofstream out(filename.c_str());
        out << "OFF\n";
        out << delaunay_->number_of_vertices() << ' '
            << delaunay_->number_of_finite_facets() << " 0\n";

        std::map<Delaunay::Vertex_handle, int> V;
        int index = 0;

        // Write the points
        for (Delaunay::Finite_vertices_iterator vit = delaunay_->finite_vertices_begin();
             vit != delaunay_->finite_vertices_end(); ++vit)
        {
            V[vit] = index++;
            Point p = vit->point();
            out << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
        }

        // Write the facets
        for (Delaunay::Finite_facets_iterator fit = delaunay_->finite_facets_begin();
             fit != delaunay_->finite_facets_end(); ++fit)
        {
            Delaunay::Cell_handle ch = fit->first;
            int i = fit->second;
            out << "3 ";
            for (int j = 0; j < 4; ++j)
            {
                if (j != i)
                {
                    out << V[ch->vertex(j)] << ' ';
                }
            }
            out << '\n';
        }
        out.close();
    }

    void MCMT::export_surface_obj(const std::string &filename)
    {

        std::vector<Point> mesh_vertices;
        std::vector<std::vector<int>> mesh_faces;
        std::map<std::pair<int, int>, int> ev_map;
        std::vector<Cell_handle> intersection_cell;

        for (Finite_cells_iterator cell_it = delaunay_->finite_cells_begin(); cell_it != delaunay_->finite_cells_end(); cell_it++)
        {
            //   std::cout << "cell_it: " << cell_it->vertex(0)->info().point_index << std::endl;

            std::vector<int> v_indices;
            bool flag = false;
            for (size_t lv = 0; lv < 4; ++lv)
            {
                int v = cell_it->vertex(lv)->info().point_index;
                v_indices.push_back(int(v));
            }

            for (std::pair vertex_pair : configurations_)
            {

                if (cell_it->vertex(vertex_pair.first)->info().point_value * cell_it->vertex(vertex_pair.second)->info().point_value < 0)
                {
                    flag = true;
                    std::pair<int, int> e01 = (v_indices[vertex_pair.first] < v_indices[vertex_pair.second]) ? std::make_pair(v_indices[vertex_pair.first], v_indices[vertex_pair.second]) : std::make_pair(v_indices[vertex_pair.second], v_indices[vertex_pair.first]);
                    if (ev_map.find(e01) == ev_map.end())
                    {
                        Point interpolated_point = interpolate(cell_it->vertex(vertex_pair.first), cell_it->vertex(vertex_pair.second));
                        ev_map.insert({e01, mesh_vertices.size()});
                        mesh_vertices.push_back(interpolated_point);
                    }
                }
            }
            if (flag)
            {
                intersection_cell.push_back(cell_it);
            }
        }
        std::cout << "Number of intersection cells " << intersection_cell.size() << std::endl;

        for (int i = 0; i < intersection_cell.size(); i++)
        {
            Cell_handle ch = intersection_cell[i];
            unsigned char index = 0;
            std::vector<int> v_indices;
            for (int lv = 0; lv < 4; ++lv)
            {
                int v = ch->vertex(lv)->info().point_index;
                // int v = delaunay_->cell_vertex(cell_index, lv);
                v_indices.push_back(int(v));

                if (ch->vertex(lv)->info().point_value < 0)
                {
                    index |= (1 << lv);
                }
            }

            std::pair<int, int> e01 = (v_indices[0] < v_indices[1]) ? std::make_pair(v_indices[0], v_indices[1]) : std::make_pair(v_indices[1], v_indices[0]);
            std::pair<int, int> e02 = (v_indices[0] < v_indices[2]) ? std::make_pair(v_indices[0], v_indices[2]) : std::make_pair(v_indices[2], v_indices[0]);
            std::pair<int, int> e03 = (v_indices[0] < v_indices[3]) ? std::make_pair(v_indices[0], v_indices[3]) : std::make_pair(v_indices[3], v_indices[0]);
            std::pair<int, int> e12 = (v_indices[1] < v_indices[2]) ? std::make_pair(v_indices[1], v_indices[2]) : std::make_pair(v_indices[2], v_indices[1]);
            std::pair<int, int> e13 = (v_indices[1] < v_indices[3]) ? std::make_pair(v_indices[1], v_indices[3]) : std::make_pair(v_indices[3], v_indices[1]);
            std::pair<int, int> e23 = (v_indices[2] < v_indices[3]) ? std::make_pair(v_indices[2], v_indices[3]) : std::make_pair(v_indices[3], v_indices[2]);

            int v00, v10, v20, v01, v11, v21;

            // Case analysis
            switch (index)
            {
                // skip inside or outside situitions
            case 0x00:
            case 0x0F:
                break;
                // only vert 0 is inside
            case 0x01:
                v00 = ev_map[e01];
                v10 = ev_map[e03];
                v20 = ev_map[e02];

                mesh_faces.push_back({v00, v20, v10});
                break;

                // only vert 1 is inside
            case 0x02:
                v00 = ev_map[e01];
                v10 = ev_map[e12];
                v20 = ev_map[e13];

                mesh_faces.push_back({v00, v20, v10});
                break;

                // only vert 2 is inside
            case 0x04:
                v00 = ev_map[e02];
                v10 = ev_map[e23];
                v20 = ev_map[e12];

                mesh_faces.push_back({v00, v20, v10});
                break;

                // only vert 3 is inside
            case 0x08:
                v00 = ev_map[e13];
                v10 = ev_map[e23];
                v20 = ev_map[e03];

                mesh_faces.push_back({v00, v20, v10});
                break;

                // verts 0, 1 are inside
            case 0x03:
                v00 = ev_map[e03];
                v10 = ev_map[e02];
                v20 = ev_map[e13];

                mesh_faces.push_back({v00, v20, v10});
                v01 = ev_map[e02];
                v11 = ev_map[e12];
                v21 = ev_map[e13];

                mesh_faces.push_back({v01, v21, v11});
                break;

                // verts 0, 2 are inside
            case 0x05:
                v00 = ev_map[e03];
                v10 = ev_map[e12];
                v20 = ev_map[e01];

                mesh_faces.push_back({v00, v20, v10});
                v01 = ev_map[e12];
                v11 = ev_map[e03];
                v21 = ev_map[e23];

                mesh_faces.push_back({v01, v21, v11});
                break;

                // verts 0, 3 are inside
            case 0x09:
                v00 = ev_map[e01];
                v10 = ev_map[e13];
                v20 = ev_map[e02];

                mesh_faces.push_back({v00, v20, v10});
                v01 = ev_map[e13];
                v11 = ev_map[e23];
                v21 = ev_map[e02];

                mesh_faces.push_back({v01, v21, v11});
                break;

                // verts 1, 2 are inside
            case 0x06:
                v00 = ev_map[e01];
                v10 = ev_map[e02];
                v20 = ev_map[e13];

                mesh_faces.push_back({v00, v20, v10});
                v01 = ev_map[e13];
                v11 = ev_map[e02];
                v21 = ev_map[e23];
                mesh_faces.push_back({v01, v21, v11});
                break;

                // verts 2, 3 are inside
            case 0x0C:
                v00 = ev_map[e13];
                v10 = ev_map[e02];
                v20 = ev_map[e03];

                mesh_faces.push_back({v00, v20, v10});
                v01 = ev_map[e02];
                v11 = ev_map[e13];
                v21 = ev_map[e12];
                mesh_faces.push_back({v01, v21, v11});
                break;

                // verts 1, 3 are inside
            case 0x0A:
                v00 = ev_map[e03];
                v10 = ev_map[e01];
                v20 = ev_map[e12];

                mesh_faces.push_back({v00, v20, v10});
                v01 = ev_map[e12];
                v11 = ev_map[e23];
                v21 = ev_map[e03];
                mesh_faces.push_back({v01, v21, v11});
                break;

                // verts 0, 1, 2 are inside
            case 0x07:
                v00 = ev_map[e03];
                v10 = ev_map[e23];
                v20 = ev_map[e13];

                mesh_faces.push_back({v00, v20, v10});
                break;

                // verts 0, 1, 3 are inside
            case 0x0B:
                v00 = ev_map[e12];
                v10 = ev_map[e23];
                v20 = ev_map[e02];

                mesh_faces.push_back({v00, v20, v10});
                break;

                // verts 0, 2, 3 are inside
            case 0x0D:
                v00 = ev_map[e01];
                v10 = ev_map[e13];
                v20 = ev_map[e12];

                mesh_faces.push_back({v00, v20, v10});
                break;

                // verts 1, 2, 3 are inside
            case 0x0E:
                v00 = ev_map[e01];
                v10 = ev_map[e02];
                v20 = ev_map[e03];

                mesh_faces.push_back({v00, v20, v10});
                break;

            default:
                // assert(false);
                continue;
            }
        }

        std::ofstream outfile(filename); // create a file named "example.txt"

        if (outfile.is_open())
        {
            for (size_t k = 0; k < mesh_vertices.size(); ++k)
            {
                Point site = mesh_vertices[k];
                // size_t n = site.size();
                outfile << "v ";
                outfile << site.x() << " " << site.y() << " " << site.z();
                outfile << " \n";
            }

            for (size_t j = 0; j < mesh_faces.size(); j++)
            {
                std::vector<int> vertices_idx = mesh_faces[j];

                outfile << "f " << vertices_idx[0] + 1 << " " << vertices_idx[1] + 1 << " " << vertices_idx[2] + 1 << " \n";
            }
        }
    }

    Point MCMT::interpolate(Vertex_handle p1, Vertex_handle p2)
    {
        double p1_x = p1->point().x();
        double p1_y = p1->point().y();
        double p1_z = p1->point().z();
        double p2_x = p2->point().x();
        double p2_y = p2->point().y();
        double p2_z = p2->point().z();
        double sd1 = p1->info().point_value;
        double sd2 = p2->info().point_value;
        double t = sd1 / ((sd1 - sd2));
        if (abs(sd1 - sd2) < 1e-6)
        {
            // std::cout << "WARNING! SD1 == SD2" << std::endl;
            t = 0.5;
        }
        return Point(p1_x + t * (p2_x - p1_x), p1_y + t * (p2_y - p1_y), p1_z + t * (p2_z - p1_z));
    }

    std::vector<double> MCMT::compute_point_bbox(const std::vector<Point> &points)
    {
        // compute min max of points
        double min_x = points[0].x();
        double min_y = points[0].y();
        double min_z = points[0].z();
        double max_x = points[0].x();
        double max_y = points[0].y();
        double max_z = points[0].z();
        for (int i = 0; i < points.size(); i++)
        {
            if (points[i].x() < min_x)
            {
                min_x = points[i].x();
            }
            if (points[i].y() < min_y)
            {
                min_y = points[i].y();
            }
            if (points[i].z() < min_z)
            {
                min_z = points[i].z();
            }
            if (points[i].x() > max_x)
            {
                max_x = points[i].x();
            }
            if (points[i].y() > max_y)
            {
                max_y = points[i].y();
            }
            if (points[i].z() > max_z)
            {
                max_z = points[i].z();
            }
        }
        double padding = 1e-3;
        min_x -= padding;
        min_y -= padding;
        min_z -= padding;
        max_x += padding;
        max_y += padding;
        max_z += padding;
        return std::vector<double>{min_x, min_y, min_z, max_x, max_y, max_z};
    }

    double MCMT::compute_point_density(double value)
    {
        return 1 / (abs(value) + 1e-6);
    }
}
