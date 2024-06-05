#pragma once

#include "Mesh.hpp"
#include "parameters.hpp"

class Solver : public Mesh{
  /**
   * @brief Class to handle the solver of the PDE
   */

  // variables for MPI to avoid re-calculation
  std::vector<int> send_counts;

  public:
  Solver(std::vector<double> & _mesh, const Domain & d, const size_t & n_col, const std::string & f);
  Solver(Mesh & m);

  void print_mesh() const;
  std::optional<std::vector<double>> solution_finder_sequential();
  void initial_communication(std::vector<double> & initial_mesh);
  void solution_finder_mpi(std::vector<double> & final_mesh, const int & thread = 4);
  void communicate_boundary();
  void final_communication(std::vector<double> & final_mesh);
};