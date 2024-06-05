#include "Solver.hpp"
#include <array>
#include <numeric>

Solver::Solver(std::vector<double> & _mesh, const Domain & d, const size_t & n_col, const std::string & f) : mesh(_mesh), mesh_obj(_mesh, n_col, d, f), n(n_col), f(f) {
  /**
   * @brief Constructor of the Solver class
   * @param _mesh is the mesh to be solved
   * @param d is the domain of the mesh
   * @param n_col is the number of points of column of the mesh
   * @param f is the function to be computed
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  send_counts.reserve(size);
}

Solver::Solver(Mesh & m, const size_t & _n) : mesh(m.get_mesh()), mesh_obj(m), n(_n), f(m.get_f()) {
  /**
   * @brief Constructor of the Solver class
   * @param m is the mesh to be solved
   * @param _n is the number of points of column of the mesh
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  send_counts.reserve(size);
}

void Solver::print_mesh() const{
  /**
   * @brief Function to print the mesh and rank who is printing
  */

  //sleep(rank); // to print in the correct order (not necessary for the final version)

  std::cout << "Rank: " << rank << std::endl;
  int spacing = 7;
  for (size_t r = 0; r < mesh.size()/n; ++r) {
    std::cout << std::setw(spacing) << r << "| ";
    for (size_t c = 0; c < n; ++c) {
        std::cout << std::setw(spacing) << std::setprecision (2) << mesh[r * n + c] << " "; // Ensure alignment for mesh values
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

std::optional<std::vector<double>> Solver::solution_finder_sequential(){
  /**
   * @brief Function to find the solution of the mesh totally sequentially
  */

  // Variables creation
  conditions c;
  bool exit = false;
  double mean_time = 0;
  int i = 1;

  do{
    // update mesh + time measure
    auto start = std::chrono::high_resolution_clock::now();
    mesh_obj.update_seq();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    
    // add the value to the mean time
    mean_time += duration.count();
    
    // check stopping criteria
    exit = mesh_obj.get_error() < c.tolerance || i == c.n_max - 1;

    if(exit)
      std::cout << "Iter: " << i<<" - time: " << mean_time << " ms - Mean time each update: "<< mean_time/i << " ms" << std::endl;

    ++i;
  }while (!exit);

  // save the final mesh into a file
  std::string filename = "vtk_files/approx_sol-1-" + std::to_string(mesh_obj.get_size().first) + ".vtk";
  mesh_obj.write(filename);

  #if TEST == 1
  return mesh_obj.get_mesh();
  #endif

  return std::nullopt;
}

void Solver::solution_finder_parallel(int n_tasks){
  /**
   * @brief Function to find the solution of the mesh using openMP
   * @param n_tasks is the number of parallel tasks - OpenMP
  */

  // Variables creation
  conditions c;
  bool exit = false;
  double mean_time = 0;
  int i = 1;

  do{
    // update mesh + time measure
    auto start = std::chrono::high_resolution_clock::now();
    mesh_obj.update_par(n_tasks);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    
    // add the value to the mean time
    mean_time += duration.count();
    
    // check stopping criteria
    exit = mesh_obj.get_error() < c.tolerance || i == c.n_max - 1;

    if(exit)
      std::cout << "Iter: " << i<<" - time: " << mean_time << " ms - Mean time each update: "<< mean_time/i << " ms"<< std::endl;

    ++i;
  }while (!exit);

  // save the final mesh into a file
  std::string filename = "vtk_files/approx_sol-1-" + std::to_string(mesh_obj.get_size().first) + ".vtk";
  mesh_obj.write(filename);
}

void Solver::communicate_boundary() {
  /**
   * @brief Function to communicate the boundary of mesh inside each MPI process
   */

  // take the last mesh
  mesh = mesh_obj.get_mesh();

  if (rank == 0) {
      // Send the last computed row to the next process and receive the first computed row from process 1
      MPI_Sendrecv(&mesh[mesh.size() - 2*n], n, MPI_DOUBLE, 1, 0, &mesh[mesh.size() - n], n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else if (rank == size - 1) {
      // Send the first computed row to the previous process and receive the last computed row from the previous process
      MPI_Sendrecv(&mesh[n], n, MPI_DOUBLE, rank - 1, 0, &mesh[0], n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
      // Send the first computed row to the previous process and receive the last computed row from the previous process
      MPI_Sendrecv(&mesh[n], n, MPI_DOUBLE, rank - 1, 0, &mesh[0], n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Send the last computed row to the next process and receive the first computed row from next process
      MPI_Sendrecv(&mesh[mesh.size() - 2*n], n, MPI_DOUBLE, rank + 1, 0, &mesh[mesh.size() - n], n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  // Save the new boundary condition inside the mesh object
  mesh_obj.set_mesh(mesh);
}

void Solver::initial_communication(std::vector<double> & initial_mesh){
  /**
   * @brief Function to communicate the initial mesh
   * @param initial_mesh is the initial mesh to be communicated
  */
  
  int row_eq_distr = (n - 2)/(size); // -2 since we have to row of border condition
  int remainder = (n - 2)%(size);

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
    displacement[i] = displacement[i-1] + temp[i-1]*n;
  }

  // offset calculation to correctly find x during calculations
  std::vector<int> offset_vec(size, 0);
  offset_vec[0] = 0; // since 0 start from 0
  std::transform(displacement.begin(), displacement.end(), offset_vec.begin(), [this](int val){ return val/n;} );
  
  /*
  if(rank == 0){
    // [DEBUG]
    for(int i = 1; i < size; ++i){offset_vec[i] = offset_vec[i-1] + temp[i-1];}
    for(auto & i: offset_vec){std::cout << i << " ";} std::endl(std::cout);
    std::transform(displacement.begin(), displacement.end(), offset_vec.begin(), [this](int val){ return val/n;} );
    for(auto & i: offset_vec){std::cout << i << " ";} std::endl(std::cout);
  }
  */

  // Offset assignment
  mesh_obj.set_offset(offset_vec[rank]);
  
  // Calculate element sent: n_row_each_process*n + 2*n which are boundaries
  send_counts.resize(size);
  std::transform(temp.begin(), temp.end(), send_counts.begin(), [this](int val){ return val != 0 ? val*n + 2*n : 0;} );
  
  // Send the initial mesh to the other threads
  MPI_Scatterv(&initial_mesh[0], &send_counts[0], &displacement[0], MPI_DOUBLE, &mesh[0], send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Save the initial condition inside the mesh object
  mesh_obj.set_mesh(mesh);
}

void Solver::final_communication(std::vector<double> & final_mesh){
  /**
   * @brief Function to communicate adn save the final mesh
   * @param mesh is the mesh to be communicated
   * @param final_mesh is the final mesh to be saved
   * @param n is the number of points of column of the mesh
   * @param f is the function to be computed
  */
  // take the last mesh
  mesh = mesh_obj.get_mesh();
  
  // Remove boundaries rows from the mesh
  std::transform(send_counts.begin(), send_counts.end(), send_counts.begin(), [this](int val){ return val != 0 ? val - 2*n : 0;} );
  
  // rows to be sent to the root process
  std::vector<int> temp(size, 0);
  std::transform(send_counts.begin(), send_counts.end(), temp.begin(), [this](int val){ return val/n;} );
  
  // calculate the displacement, offset of which root process has to write data
  std::vector<int> displacement(size, 0);
  for(int i = 1; i < size; ++i){
    displacement[i] = displacement[i-1] + send_counts[i-1];
  }

  // Starting mesh from n to avoid the first row of boundary condition
  MPI_Gatherv(&mesh[n], send_counts[rank], MPI_DOUBLE, &final_mesh[n], &send_counts[0], &displacement[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(rank == 0){
    // save the final mesh
    Mesh last_mesh(final_mesh, n, Domain(0, 1, 0, 1), f);
    //last_mesh.print();
  
    // write the final mesh
    std::string filename = "vtk_files/approx_sol-"+std::to_string(size)+"-" + std::to_string(last_mesh.get_size().first) + ".vtk";
    last_mesh.write(filename);
  }
}

void Solver::solution_finder_mpi(std::vector<double> & final_mesh, const int & thread){
  /**
   * @brief Function to find the solution of the problem using MPI + OpenMP
   * @param final_mesh is the final mesh to be communicated
   * @param thread is the number of parallel task
  */

  // Variables creation
  conditions c;
  int exit_local = mesh_obj.get_size().first == 0 ? 1 : 0; // if the mesh is empty the local exit condition is reached
  int exit = 0;
  double mean_time = 0;
  int i = 1;

  do{
    // update mesh if the local exit conditions is not reached
    if(!exit_local){

      // update mesh in parallel + measure time
      auto start = std::chrono::high_resolution_clock::now();
      mesh_obj.update_par(thread);
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

      // add the value to the mean time
      mean_time += duration.count();

      // check stopping criteria
      exit_local = (mesh_obj.get_error() < c.tolerance || i == c.n_max - 1) ? 1 : 0;
    }

    // Communicate local exit condition to all threads
    MPI_Allreduce(&exit_local, &exit, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

    // Communicate new boundary of each mesh to the other processes 
    communicate_boundary();

    // check if the global exit condition isn't reached => each process has reached its local exit condition
    if(!exit){
      // update mesh with the new boundary
      mesh_obj.set_mesh(mesh);
      
      // iter counter update
      ++i;
    }

  }while (!exit);

  // calculate the mean time
  MPI_Allreduce(MPI_IN_PLACE, &mean_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // print mean time
  if(rank == 0){
    mean_time /= size;
    std::cout << "Iter: " << i<<" - time: " << mean_time << " ms - Mean time each update: "<< mean_time/i << " mus"<< std::endl;
  }

  // Communicate result
  final_communication(final_mesh);
}
