#pragma once

#include "Mesh.hpp"

struct conditions{
  double tolerance = 1e-5;
  int n_max = 3e3;
};

class Solver{
    /**
     * @brief Class to handle the solver of the PDE
     */

    std::vector<double> mesh;
    Mesh mesh_obj;
    size_t n;
    std::string f;
    int rank, size;

    public:
    Solver(std::vector<double> & _mesh, const Domain & d, const size_t & n_col, const std::string & f);
    Solver(Mesh & m, const size_t & _n);

    void print_mesh() const;
    void solution_finder_sequential();
    void solution_finder_parallel(int n_tasks = 4);
    void communicate_boundary();
    void initial_communication(std::vector<double> & initial_mesh);
    void final_communication(std::vector<double> & final_mesh);
    void solution_finder_mpi(std::vector<double> & final_mesh, const int & thread = 4);

};