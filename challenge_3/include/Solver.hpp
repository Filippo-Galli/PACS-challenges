#pragma once

#include "Mesh.hpp"
#include<iostream>
#include<cmath>
#include<vector>
#include <mpi.h>

struct conditions{
  double tolerance = 1e-6;
  int n_max = 3e3;
};

class Solver{
    /**
     * @brief Class to solve the PDE
     */

    public:
    void print_mesh(const std::vector<double> & mesh, const size_t & n);
    void solution_finder_sequential(Mesh & mesh, int n_tasks = 4);
    void communicate_boundary(std::vector<double> & mesh, const size_t & n);
    void initial_communication(std::vector<double> & initial_mesh, const size_t & n, std::vector<double> & mesh);
    void final_communication(std::vector<double> & mesh, std::vector<double> & final_mesh, const int & n, const std::string & f);
    void solution_finder_mpi(std::vector<double> & mesh, const size_t & n, const std::string & f, std::vector<double> & final_mesh, Domain d, const int & thread = 4);

};