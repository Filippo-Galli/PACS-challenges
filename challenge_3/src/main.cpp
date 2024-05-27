
#include<iostream>
#include<cmath>
#include<vector>

#include "Solver.hpp"

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

    Solver sol;

    sol.solution_finder_sequential(mesh, atoi(argv[3]));

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

    Solver sol;

    // communicate the initial mesh
    sol.initial_communication(temp_rank0, n, mesh);

    sol.solution_finder_mpi(mesh, atoi(argv[1]), argv[2], temp_rank0, domain, atoi(argv[3]));
  }

  MPI_Finalize();
  return 0;
}
