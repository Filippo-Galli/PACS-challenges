#include "Solver.hpp"

bool check_input(int argc, char *argv[]){
  // check the number of arguments
  if(argc != 5) {
    std::cout << "Usage: " << argv[0] << " [point of the mesh]" << " [function]" << " [number of parallel task]"  << " [function of boundaries]"<< std::endl;
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

void distance_from_exact_solution(const Mesh & mesh){
  /**
   * @brief Function to calculate the distance from the exact solution to test this program
   * @param mesh is the mesh object of the approximate solution
   * 
  */

  conditions cond;
  
  double error = 0;
  
  for(size_t i = 0; i < mesh.get_size().first; ++i){
    for(size_t j = 0; j < mesh.get_size().second; ++j){
      auto [x, y] = mesh.get_coordinates(i, j);
      error += pow(mesh.get_value(i, j) - cond.u(x, y), 2);
    }
  }

  std::cout << "Error with exact solution: " << sqrt(mesh.get_h()*error) << std::endl;
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
  
  size_t n = atoi(argv[1]);
  Domain domain(0, 1, 0, 1);

  if(size == 1){
    // create a mesh
    Mesh mesh(n, n, domain, argv[2]);

    // cdd boundary condition
    mesh.add_boundary_condition(argv[4]);

    // create the solver object
    Solver sol(mesh);

    // find the solution
    auto solution_vec = sol.solution_finder_sequential();

    #if TEST == 1
    // run it only if build with test flag
    distance_from_exact_solution(Mesh(solution_vec.value(), n, domain, argv[2]));
    #endif

  }
  else{
    
    // calculate the number of rows for each thread
    int row_eq_distr = (n - 2)/(size); // -2 since we have to row of border condition
    int remainder = (n - 2)%(size);

    std::vector<int> temp(size, row_eq_distr);
    
    // distribute the remainder
    int idx = 0;
    for(; remainder > 0; --remainder, ++idx){
      ++temp[idx];
      if(idx == size)
        idx = 0; 
    }

    // declare variable
    std::vector<double> mesh_vec(temp[rank]*n != 0? temp[rank]*n + 2*n : 3*n, 0);
    std::vector<double> total_mesh;

    // correctly resize meshes
    if(rank == 0){
      total_mesh = std::vector<double>(n*n, 0); 

      Mesh templ_mesh(total_mesh, n, domain, argv[2]);

      // add border condition if needed
      templ_mesh.add_boundary_condition(argv[4]);

      // copy the mesh to the total mesh with border conditions
      total_mesh.assign(templ_mesh.get_mesh().begin(), templ_mesh.get_mesh().end());
    }

    // initialize solver
    Solver sol(mesh_vec, domain, n, argv[2]);

    // initial communication
    sol.initial_communication(total_mesh);

    // find the solution
    sol.solution_finder_mpi(total_mesh, atoi(argv[3]));    

    #if TEST == 1
    // run it only if build with test flag
    if(rank == 0){
      distance_from_exact_solution(Mesh(total_mesh, n, domain, argv[2]));
    }
    #endif
  }
  
  MPI_Finalize();
  return 0;
}
