#include "Solver.hpp"

void Solver::print_mesh() const{
  /**
   * @brief Function to print the mesh and rank who is printing
   * @param mesh is the mesh to be printed
   * @param n is the number of points of column of the mesh
  */
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

void Solver::solution_finder_sequential(){
  /**
   * @brief Function to find the solution of the mesh
   * @param mesh is the mesh to be solved
  */

  // update the mesh until convergence
  conditions c;
  bool exit = false;

  double mean_time = 0;

  int i = 1;
  do{
    // update mesh
    auto start = std::chrono::high_resolution_clock::now();

    mesh_obj.update_seq();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    // add the value to the mean time
    mean_time += duration.count();
    
    // check stopping criteria
    exit = mesh_obj.get_error() < c.tolerance || i == c.n_max - 1;

    if(exit)
      std::cout << "Iter: " << i<<" - time: " << mean_time/1e6 << " s - Mean time each update: "<< mean_time/i << " mus - error: " << std::setprecision(9) << mesh_obj.get_error() << std::endl;

    ++i;
  }while (!exit);

  // save the final mesh
  std::string filename = "vtk_files/approx_sol-1-" + std::to_string(mesh_obj.get_size().first) + ".vtk";
  mesh_obj.write(filename);
}

void Solver::solution_finder_parallel(int n_tasks){
  /**
   * @brief Function to find the solution of the mesh
   * @param mesh is the mesh to be solved
  */

  // update the mesh until convergence
  conditions c;
  bool exit = false;

  double mean_time = 0;

  int i = 1;
  do{
    // update mesh
    auto start = std::chrono::high_resolution_clock::now();

    mesh_obj.update_par(n_tasks);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    // add the value to the mean time
    mean_time += duration.count();
    
    // check stopping criteria
    exit = mesh_obj.get_error() < c.tolerance || i == c.n_max - 1;

    if(exit)
      std::cout << "Iter: " << i<<" - time: " << mean_time/1e6 << " s - Mean time each update: "<< mean_time/i << " mus - error: " << std::setprecision(9) << mesh_obj.get_error() << std::endl;

    ++i;
  }while (!exit);

  // save the final mesh
  std::string filename = "vtk_files/approx_sol-1-" + std::to_string(mesh_obj.get_size().first) + ".vtk";
  mesh_obj.write(filename);
}

void Solver::communicate_boundary(){
  /**
   * @brief Function to communicate the boundary of the mesh
   * @param n is the number of points of column of the mesh
  */

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank == 0){
    // send the last computed row to the next thread
    MPI_Send(&mesh[mesh.size() - 2*n], n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    // receive the last computed row from thread 1
    MPI_Recv(&mesh[mesh.size() - n], n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } 
  else if(rank == size - 1){
    MPI_Send(&mesh[n], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD);
    MPI_Recv(&mesh[0], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else{
    MPI_Send(&mesh[n], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD);
    MPI_Recv(&mesh[0], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    MPI_Recv(&mesh[mesh.size() - n], n, MPI_DOUBLE, rank +1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&mesh[mesh.size() - 2*n], n, MPI_DOUBLE, rank +1, 0, MPI_COMM_WORLD);
  }
}

void Solver::initial_communication(std::vector<double> & initial_mesh){
  /**
   * @brief Function to communicate the initial mesh
   * @param initial_mesh is the initial mesh to be communicated
   * @param n is the number of points of column of the mesh
   * @param mesh is the mesh of each thread, the receiving buffer
  */

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t expected_size;
  if(rank == 0 or rank == size - 1)
    expected_size = n*n/size + n;
  else
    expected_size = n*n/size + 2*n;

  if(rank == 0){
    // select only first 3 row of initial mesh to save into mesh 
    for(size_t i = 0; i < expected_size; ++i)
      mesh[i] = initial_mesh[i];
    
    for(int proc = 1; proc < size; ++proc){
      if(proc == size - 1)
        MPI_Send(&initial_mesh[(proc*n/size -1)*n], n*n/size + n, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
      else
        MPI_Send(&initial_mesh[(proc*n/size -1)*n], expected_size + n, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(&mesh[0], expected_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

}

void Solver::final_communication(std::vector<double> & final_mesh){
  /**
   * @brief Function to communicate adn save the final mesh
   * @param mesh is the mesh to be communicated
   * @param final_mesh is the final mesh to be saved
   * @param n is the number of points of column of the mesh
   * @param f is the function to be computed
  */

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if(rank == 0){
    MPI_Gather(&mesh[0], mesh.size() - n, MPI_DOUBLE, &final_mesh[0], mesh.size()  - n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else if (rank == size -1){
    MPI_Gather(&mesh[n], mesh.size() - n, MPI_DOUBLE, &final_mesh[0], mesh.size()  - n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else{
    mesh.resize(mesh.size() - n); // delete last row
    MPI_Gather(&mesh[n], mesh.size() -n, MPI_DOUBLE, &final_mesh[0], mesh.size()  - n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  if(rank == 0){
    // save the final mesh
    Mesh last_mesh(final_mesh, n, Domain(0, 1, 0, 1), f);
    
    // write the final mesh
    std::string filename = "vtk_files/approx_sol-"+std::to_string(size)+"-" + std::to_string(last_mesh.get_size().first) + ".vtk";
    last_mesh.write(filename);
  }
}

void Solver::solution_finder_mpi(std::vector<double> & final_mesh, const int & thread){
  /**
   * @brief Function to find the solution of the problem
   * @param mesh is the mesh to be computed
   * @param n is the number of points of column of the mesh
   * @param f is the function to be computed
   * @param final_mesh is the final mesh to be communicated
   * @param thread is the number of parallel task
  */

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  conditions c;
  int exit = 0;
  int exit_local = 0;

  double mean_time = 0;

  int i = 1;
  do{
    // update mesh if the exit conditions is not reached
    if(!exit_local){

      // update mesh in parallel + measure time
      auto start = std::chrono::high_resolution_clock::now();
      mesh_obj.update_par(thread);
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

      mean_time += duration.count();

      // check stopping criteria
      exit_local = (mesh_obj.get_error() < c.tolerance || i == c.n_max - 1) ? 1 : 0;
    }

    // Communicate exit condition to all threads
    MPI_Allreduce(&exit_local, &exit, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

    mesh = mesh_obj.get_mesh();
    // Communicate mesh 
    communicate_boundary();
    
    if(exit){
      std::cout << "Rank: "<< rank << " - iter: " << i<<" - time: " << std::setw(8) << mean_time/1e6 << " s - Mean time each update: "<< std::setw(8)<< mean_time/i << " mus - error: " << std::setprecision(9) << mesh_obj.get_error() << std::endl;
    }
    else{
      // update mesh with the new boundary
      mesh_obj.set_mesh(mesh);
      ++i;
    }

  }while (!exit);

  // Communicate result
  mesh = mesh_obj.get_mesh();
  final_communication(final_mesh);
  
}