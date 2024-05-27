#include "Solver.hpp"

bool check_input(int argc, char *argv[]){
  // check the number of arguments
  if(argc != 4) {
    std::cout << "Usage: " << argv[0] << " [point of the mesh]" << " [function]" << " [number of parallel task]" << std::endl;
    return false;
  }

  // check on the number of points of the mesh
  if(atoi(argv[1]) < 1) {
    std::cout << "The number of points must be greater than 0" << std::endl;
    return false;
  }
  
  // check on the number of parallel tasks
  if(atoi(argv[3]) < 1) {
    std::cout << "The number of parallel tasks must be greater than 0" << std::endl;
    return false;
  }

  return true;
}


int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(!check_input(argc, argv)){
    MPI_Finalize();
    return 1;
  }
  
  int n = atoi(argv[1]);
  Domain domain(0, 1, 0, 1);

  if(size == 1){
    // create a mesh
    Mesh mesh(n, n, domain, argv[2]);

    // Create the solver object
    Solver sol(mesh, n);

    // find the solution
    sol.solution_finder_sequential();

  }
  else{
    
    // declare sub-mesh and total one
    std::vector<double> mesh_vec(n*n/size);
    std::vector<double> total_mesh;

    // correctly resize meshes
    if(rank == 0){
      total_mesh = std::vector<double>(n*n, 0);
      mesh_vec.resize(n*n/size + n); // +1 row since it has only 1 row as "boundary"
    }
    else if(rank == size - 1)
      mesh_vec.resize(n*n/size + n); // +1 row since it has only 1 row as "boundary"
    else
      mesh_vec.resize(n*n/size + 2*n); // +2 row which are the ones of the other threads

    // Initialize solver
    Solver sol(mesh_vec, domain, n, argv[2]);

    // communicate total mesh
    sol.initial_communication(total_mesh);

    // find the solution
    sol.solution_finder_mpi(total_mesh, atoi(argv[3]));
  }

  MPI_Finalize();
  return 0;
}
