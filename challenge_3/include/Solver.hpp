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

    public:
    Solver(std::vector<double> & _mesh, const Domain & d, const size_t & n_col, const std::string & f) : mesh(_mesh), mesh_obj(_mesh, n_col, d, f), n(n_col), f(f) {};
    Solver(Mesh & m, const size_t & _n) : mesh(m.get_mesh()), mesh_obj(m), n(_n), f(m.get_f()) {};

    void print_mesh() const;
    void solution_finder_sequential();
    void solution_finder_parallel(int n_tasks = 4);
    void communicate_boundary();
    void initial_communication(std::vector<double> & initial_mesh);
    void final_communication(std::vector<double> & final_mesh);
    void solution_finder_mpi(std::vector<double> & final_mesh, const int & thread = 4);

};