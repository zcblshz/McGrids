#include <chrono>
#include <fstream>
#include <sstream>
#include <set>
#include <iostream>
#include <vector>
#include <tbb/tbb.h>
#include "KDTreeVectorOfVectorsAdaptor.hpp"
#define MEASUREMENTS 600
#define CUT_MEASUREMENTS 100
#define GB 1000000000

class KDTree
{
public:
    KDTree(int num_points, double *point_positions, double *point_values)
    {
        for (int i = 0; i < num_points; i++)
        {
            kdtree_point_values.push_back(point_values[i]);
            std::vector<double> point;
            for (int j = 0; j < 3; j++)
            {
                point.push_back(point_positions[i * 3 + j]);
            }
            point.push_back(point_values[i]);
            kdtree_point_positions.push_back(point);
        }
    }
    ~KDTree()
    {
        kdtree_point_positions.clear();
        kdtree_point_values.clear();
    }

    std::vector<double> compute_density(int num_points, double *point_positions)
    {
        typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<double>>, double> my_kd_tree_t;

        my_kd_tree_t tree(3, kdtree_point_positions, 10 /* max leaf */);
        tree.index->buildIndex();
        std::vector<double> densities;
        for (int i = 0; i < num_points; i++)
        {
            nanoflann::KNNResultSet<double> resultSet(1);
            std::vector<size_t> ret_indexes(1);
            std::vector<double> out_dists_sqr(1);
            resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
            tree.index->findNeighbors(resultSet, point_positions + i * 3, nanoflann::SearchParameters(10));
            densities.push_back(kdtree_point_values[ret_indexes[0]]);
        }
        return densities;
    }

private:
    std::vector<std::vector<double>> kdtree_point_positions;
    std::vector<double> kdtree_point_values;
    int num_points;
};
