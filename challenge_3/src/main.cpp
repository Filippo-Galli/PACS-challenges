#include"Mesh.hpp"

void solution_finder(Mesh & mesh, int n_tasks = 4){
  /**
   * @brief Function to find the solution of the mesh
   * @param mesh is the mesh to be solved
  */

  // update the mesh until convergence
  double tolerance = 1e-6;
  int n_max = 3e3;
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
    
    exit = mesh.get_error() < tolerance || i == n_max - 1;
    if(exit)
      std::cout << "Iter: " << i<<" - time: " << mean_time/1e6 << " s - Mean time each update: "<< mean_time/i << " mus - error: " << std::setprecision(9) << mesh.get_error() << std::endl;

    ++i;
  }while (!exit);
}

int main(int argc, char *argv[]) {

  /*
  MPI_Init(&argc, &argv);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
          processor_name, world_rank, world_size);
  */
  
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
  
  Domain domain(0, 1, 0, 1);

  // create a mesh
  Mesh mesh(atoi(argv[1]), domain, argv[2]);

  solution_finder(mesh, atoi(argv[3]));

  mesh.write("vtk_files/approx_sol.vtk");
  

  //MPI_Finalize();

  return 0;
}
