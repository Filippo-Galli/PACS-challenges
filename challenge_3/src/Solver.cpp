#include "Solver.hpp"
#include <array>
#include <numeric>

Solver::Solver(std::vector<double> & _mesh, const Domain & d, const size_t & n_col, const std::string & f) : Mesh(_mesh, n_col, d, f){
  /**
   * @brief Constructor of the Solver class
   * @param _mesh is the mesh to be solved
   * @param d is the domain of the mesh
   * @param n_col is the number of points of column of the mesh
   * @param f is the function to be computed
  */

  send_counts.reserve(size);
}

Solver::Solver(Mesh & m) : Mesh(m) {
  /**
   * @brief Constructor of the Solver class
   * @param m is the mesh to be solved
   * @param _n is the number of points of column of the mesh
  */

  send_counts.reserve(size);
}

void Solver::print_mesh() const{
  /**
   * @brief Function to print the mesh and rank who is printing
  */

  std::cout << "Rank: " << rank << std::endl;
  int spacing = 7;
  for (size_t r = 0; r < mesh.size()/n_col; ++r) {
    std::cout << std::setw(spacing) << r << "| ";
    for (size_t c = 0; c < n_col; ++c) {
        std::cout << std::setw(spacing) << std::setprecision (2) << mesh[r * n_col + c] << " "; // Ensure alignment for mesh values
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

std::optional<std::vector<double>> Solver::solution_finder_sequential(){
  /**
   * @brief Function to find the solution of the mesh totally sequentially
  */

  // Create evaluation of f
  f_eval_creation();

  // Variables creation
  conditions c;
  int iter = 0;
  double e = 10;

  // Sequential computation
  auto start = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < c.n_max && e > c.tolerance; ++i){
    update_par();
    ++iter;
    e = get_error();
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time: " << duration.count() << " ms" << " - Iter: " << iter << std::endl;

  // save the final mesh into the mesh of solver
  mesh = get_mesh();

  // save the final mesh into a file
  std::string filename = "vtk_files/approx_sol-1-" + std::to_string(n_row) + ".vtk";
  write(filename);

  #if TEST == 1
  return mesh;
  #endif

  return std::nullopt;
}

void Solver::communicate_boundary() {
  /**
   * @brief Function to communicate the boundary of mesh inside each MPI process
   */

  if (rank == 0) {
      // Send the last computed row to the next process and receive the first computed row from process 1
      MPI_Sendrecv(&mesh[mesh.size() - 2*n_col], n_col, MPI_DOUBLE, 1, 0, &mesh[mesh.size() - n_col], n_col, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else if (rank == size - 1) {
      // Send the first computed row to the previous process and receive the last computed row from the previous process
      MPI_Sendrecv(&mesh[n_col], n_col, MPI_DOUBLE, rank - 1, 0, &mesh[0], n_col, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
      // Send the first computed row to the previous process and receive the last computed row from the previous process
      MPI_Sendrecv(&mesh[n_col], n_col, MPI_DOUBLE, rank - 1, 0, &mesh[0], n_col, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Send the last computed row to the next process and receive the first computed row from next process
      MPI_Sendrecv(&mesh[mesh.size() - 2*n_col], n_col, MPI_DOUBLE, rank + 1, 0, &mesh[mesh.size() - n_col], n_col, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void Solver::initial_communication(std::vector<double> & initial_mesh){
  /**
   * @brief Function to communicate the initial mesh
   * @param initial_mesh is the initial mesh to be communicated
  */
  
  int row_eq_distr = (n_col - 2)/(size); // -2 since we have to row of border condition
  int remainder = (n_col - 2)%(size);

  // calculate the number of rows for each process equally distrubuted 
  std::vector<int> temp(size, row_eq_distr);
  
  // Distribute the remainder among processes
  int idx = 0;
  for(; remainder > 0; --remainder, ++idx){
    ++temp[idx];
    if(idx == size)
      idx = 0; 
  }

  // calculate the displacement, offset of which each process read data from
  std::vector<int> displacement(size, 0);
  for(int i = 1; i < size; ++i){
    displacement[i] = displacement[i-1] + temp[i-1]*n_col;
  }

  // offset calculation to correctly find x during calculations
  std::vector<int> offset_vec(size, 0);
  offset_vec[0] = 0; // since 0 start from 0
  std::transform(displacement.begin(), displacement.end(), offset_vec.begin(), [this](int val){ return val/n_col;} );

  // Offset assignment
  offset = offset_vec[rank];
  
  // Calculate element sent: n_row_each_process*n + 2*n which are boundaries
  send_counts.resize(size);
  std::transform(temp.begin(), temp.end(), send_counts.begin(), [this](int val){ return val != 0 ? val*n_col + 2*n_col : 0;} );
  
  // Send the initial mesh to the other threads
  MPI_Scatterv(&initial_mesh[0], &send_counts[0], &displacement[0], MPI_DOUBLE, &mesh[0], send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Solver::final_communication(std::vector<double> & final_mesh){
  /**
   * @brief Function to communicate adn save the final mesh
   * @param mesh is the mesh to be communicated
   * @param final_mesh is the final mesh to be saved
   * @param n is the number of points of column of the mesh
   * @param f is the function to be computed
  */
  
  // Remove boundaries rows from the mesh
  std::transform(send_counts.begin(), send_counts.end(), send_counts.begin(), [this](int val){ return val != 0 ? val - 2*n_col : 0;} );
  
  // rows to be sent to the root process
  std::vector<int> temp(size, 0);
  std::transform(send_counts.begin(), send_counts.end(), temp.begin(), [this](int val){ return val/n_col;} );
  
  // calculate the displacement, offset of which root process has to write data
  std::vector<int> displacement(size, 0);
  for(int i = 1; i < size; ++i){
    displacement[i] = displacement[i-1] + send_counts[i-1];
  }

  // Starting mesh from n to avoid the first row of boundary condition
  MPI_Gatherv(&mesh[n_col], send_counts[rank], MPI_DOUBLE, &final_mesh[n_col], &send_counts[0], &displacement[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(rank == 0){
    // save the final mesh
    mesh = final_mesh;
    n_row = n_col;
  
    // write the final mesh
    std::string filename = "vtk_files/approx_sol-"+std::to_string(size)+"-" + std::to_string(n_col) + ".vtk";
    write(filename);
  }
}

void Solver::solution_finder_mpi(std::vector<double> & final_mesh, const int & thread){
  /**
   * @brief Function to find the solution of the mesh using MPI
   * @param final_mesh is the final mesh to be saved
   * @param thread is the number of threads to be used by openMP
  */

  // Create evaluation of f
  f_eval_creation();

  // Variables creation
  conditions c;
  int iter = 0;
  double e = 10;
  int exit = n_row == 0 ? 1 : 0; // if the mesh is empty the local exit condition is reached

  // Parallel computation
  auto start = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < c.n_max && exit < size ; ++i){
    update_par(thread);
    ++iter;
    e = get_error();
    exit = (e < c.tolerance || i == c.n_max - 1 ) ? 1 : 0;

    // Communicate local exit condition to all threads
    MPI_Allreduce(MPI_IN_PLACE, &exit, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Communicate new boundary of each mesh to the other processes 
    communicate_boundary();
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  if(rank == 0)
    std::cout << "Time: " << duration.count() << " ms" << " - Iter: " << iter << std::endl;

  final_communication(final_mesh);
}
