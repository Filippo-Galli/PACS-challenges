#include"Mesh.hpp"

void solution_finder(Mesh & mesh){
  /**
   * @brief Function to find the solution of the mesh
   * @param mesh is the mesh to be solved
  */

  // update the mesh until convergence
  double tolerance = 1e-6;
  int n_max = 1e3;
  bool exit = false;

  double mean_time = 0;

  int i = 1;
  do{
    // update mesh
    auto start = std::chrono::high_resolution_clock::now();
    mesh.update_seq();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    mean_time += duration.count();
    
    exit = mesh.get_error() < tolerance || i == n_max - 1;
    if(exit)
      std::cout << "Iter: " << i<<" - time: " << mean_time/1e6 << " s - Mean time 1 update: "<< mean_time/i << " mus - error: " << std::setprecision(9) << mesh.get_error() << std::endl;

    ++i;
  }while (!exit);
}

int main(int argc, char *argv[]) {
  
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

  solution_finder(mesh);

  mesh.write("vtk_files/approx_sol.vtk");

  return 0;
}
