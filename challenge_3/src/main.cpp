#include"Mesh.hpp"
#include<iostream>
#include<cmath>
#include<vector>

struct conditions{
  double tolerance = 1e-6;
  int n_max = 3e3;
};

void solution_finder(Mesh & mesh, int n_tasks = 4){
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

    // mesh.update_seq();
    mesh.update_par(n_tasks);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    mean_time += duration.count();
    
    exit = mesh.get_error() < c.tolerance || i == c.n_max - 1;
    if(exit)
      std::cout << "Iter: " << i<<" - time: " << mean_time/1e6 << " s - Mean time each update: "<< mean_time/i << " mus - error: " << std::setprecision(9) << mesh.get_error() << std::endl;

    ++i;
  }while (!exit);

  mesh.write("vtk_files/approx_sol.vtk");
}

void print_mesh(const std::vector<double> & mesh, const size_t & n){
  /**
   * @brief Function to print the mesh
   * @param mesh is the mesh to be printed
   * @param n is the number of points of column of the mesh
  */
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::cout << "Rank: " << rank << std::endl;
  int spacing = 7;
  for (size_t r = 0; r < mesh.size()/n; ++r) {
    std::cout << std::setw(spacing) << r << "| "; // Use setw here as well
    for (size_t c = 0; c < n; ++c) {
        std::cout << std::setw(spacing) << std::setprecision (2) << mesh[r * n + c] << " "; // Ensure alignment for mesh values
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void communicate_boundary(std::vector<double> & mesh, const size_t & n){
  /**
   * @brief Function to communicate the boundary of the mesh
   * @param mesh is the mesh to be communicated
   * @param n is the number of points of column of the mesh
  */

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank == 0){
    MPI_Send(&mesh[mesh.size() -1 - 2*n], n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD); // send the last computed row to the next thread
    MPI_Recv(&mesh[mesh.size() - n], n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the last computed row from the next thread
  } 
  else if(rank == size - 1){
    MPI_Send(&mesh[n], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD);
    MPI_Recv(&mesh[0], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else{
    MPI_Send(&mesh[n], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD);
    MPI_Recv(&mesh[0], n, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    MPI_Recv(&mesh[mesh.size() - n], n, MPI_DOUBLE, rank +1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&mesh[mesh.size() -1 - 2*n], n, MPI_DOUBLE, rank +1, 0, MPI_COMM_WORLD);
  }
}

void initial_communication(std::vector<double> & initial_mesh, const size_t & n, std::vector<double> & mesh){
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

void solution_finder_mpi(std::vector<double> & mesh, const size_t & n, const std::string & f, std::vector<double> & final_mesh, Domain d, const int & thread = 4){
  /**
   * @brief Function to find the solution of the problem
   * @param mesh is the mesh to be computed
   * @param n is the number of points of column of the mesh
   * @param f is the function to be computed
   * @param final_mesh is the final mesh to be communicated
   * @param thread is the number of parallel task
  */

  // TODO: improve it because doesn't reach the correct result but is very close to it

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Domain domain(d.x0, d.x1, d.y0, d.y1);

  Mesh mesh_obj(mesh, n, domain, f);
  //std::cout << "[DEBUG - solution finder mpi 2] Rank: " << rank << " mesh size: " << mesh.size()<< std::endl;

  conditions c;
  int exit = 0;

  double mean_time = 0;

  int i = 1;
  do{
    // update mesh
    auto start = std::chrono::high_resolution_clock::now();
    //std::cout << "[DEBUG - solution finder mpi 2] Rank: " << rank << " mesh size: " << mesh.size()<< std::endl;
    mesh_obj.update_par(thread);
    //std::cout << "[DEBUG - solution finder mpi 3] Rank: " << rank << " mesh size: " << mesh.size()<< std::endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    mean_time += duration.count();

    exit = (mesh_obj.get_error() < c.tolerance || i == c.n_max - 1) ? 1 : 0;

    MPI_Allreduce(MPI_IN_PLACE, &exit, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    mesh = mesh_obj.get_mesh();
    communicate_boundary(mesh, n);

    if(exit == size){
      std::cout << "Iter: " << i<<" - time: " << mean_time/1e6 << " s - Mean time each update: "<< mean_time/i << " mus - error: " << std::setprecision(9) << mesh_obj.get_error() << std::endl;
    }
    else{
      mesh_obj.set_mesh(mesh);
      ++i;
    }

  }while (!exit);

  

  // Communicate result
  mesh = mesh_obj.get_mesh();
  
  if(rank == 0){
    //std::cout << "[DEBUG] Rank: " << rank << " mesh size: " << mesh.size() -n<< std::endl;
    MPI_Gather(&mesh[0], mesh.size() - n, MPI_DOUBLE, &final_mesh[0], mesh.size()  - n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else if (rank == size -1){
    //std::cout << "[DEBUG] Rank: " << rank << " mesh size: " << mesh.size() -n<< std::endl;
    MPI_Gather(&mesh[n], mesh.size() - n, MPI_DOUBLE, &final_mesh[0], mesh.size()  - n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else{
    mesh.resize(mesh.size() - n); // delete last row
    //std::cout << "[DEBUG] Rank: " << rank << " mesh size: " << mesh.size() -n<< std::endl;
    MPI_Gather(&mesh[n], mesh.size() -n, MPI_DOUBLE, &final_mesh[0], mesh.size()  - n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  // save the final mesh
  if(rank == 0){
    print_mesh(final_mesh, n);
    Mesh last_mesh(final_mesh, n, Domain(0, 1, 0, 1), f);
    
    last_mesh.write("vtk_files/approx_sol.vtk");
  }
  
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // check the number of arguments
  if(argc != 4) {
    std::cout << "Usage: " << argv[0] << " [point of the mesh]" << " [function]" << " [number of parallel task]" << std::endl;
    return 1;
  }

  // check on the number of points of the mesh
  if(atoi(argv[1]) < 1) {
    std::cout << "The number of points must be greater than 0" << std::endl;
    return 1;
  }
  
  // check on the number of parallel tasks
  if(atoi(argv[3]) < 1) {
    std::cout << "The number of parallel tasks must be greater than 0" << std::endl;
    return 1;
  }
  
  int n = atoi(argv[1]);
  Domain domain(0, 1, 0, 1);

  if(size == 1){
    // create a mesh
    Mesh mesh(n, n, domain, argv[2]);

    solution_finder(mesh, atoi(argv[3]));

    print_mesh(mesh.get_mesh(), n);

  }
  else{
    std::vector<double> mesh(n*n/size);
    std::vector<double> temp_rank0;
    if(rank == 0 or rank == size - 1)
      mesh.resize(n*n/size + n); // +1 row since it has only 1 row as "boundary"
    else
      mesh.resize(n*n/size + 2*n); // +2 row which are the ones of the other threads
    
    // create a mesh
    if(rank == 0){
      temp_rank0 = std::vector<double>(n*n, 0);      
    }

    // communicate the initial mesh
    initial_communication(temp_rank0, n, mesh);

    solution_finder_mpi(mesh, atoi(argv[1]), argv[2], temp_rank0, domain, atoi(argv[3]));
  }
  
  /*
  auto u = std::function<double(double, double)>([](double x, double y) -> double {return sin(2*std::numbers::pi*x)*sin(2*std::numbers::pi*y);});
  double error = 0;

  for(size_t i = 0; i < mesh.size(); ++i)
    for(size_t j = 0; j < mesh.size(); ++j)
      error += std::pow(mesh.get(i, j).value() - u(mesh.get_coordinates(i, j).first, mesh.get_coordinates(i, j).second), 2);

  error = std::sqrt(mesh.get_mesh_spacing()* error);
  
  std::cout << "Error with exact solution is: " << error << std::endl;
  */

  MPI_Finalize();
  return 0;
}
