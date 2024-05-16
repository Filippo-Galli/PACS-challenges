#include"Mesh.hpp"

void dist_exact_sol(const Mesh & mesh, auto u) {
  /**
   * @brief Function to compute the error between the mesh and the exact solution
   * @param mesh is the mesh to be compared
   * @param u is the exact solution
  */
  double error = 0;
  for(size_t i = 0; i < mesh.size(); ++i) {
    for(size_t j = 0; j < mesh.size(); ++j) {
      auto coords = mesh.get_coordinates(i, j);
      error += (mesh.get(i, j).value() - u(coords.first, coords.second))*(mesh.get(i, j).value() - u(coords.first, coords.second));
    }
  }
  std::cout << "Error between mesh and correct solution: " << std::setprecision(9)<< sqrt(error) << std::endl;
}

void solution_finder(Mesh & mesh, auto f){
  /**
   * @brief Function to find the solution of the mesh
   * @param mesh is the mesh to be solved
  */

  // update the mesh until convergence
  double tolerance = 1e-6;
  int n_max = 1e5;
  bool exit = false;

  int i = 0;
  do{
    // update mesh
    auto start = std::chrono::high_resolution_clock::now();
    mesh.update_seq(f);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    exit = mesh.get_error() < tolerance || i == n_max - 1;
    if(exit)
      std::cout << "Iter: " << i<<" - time: " << duration.count() << " mus - error: " << std::setprecision(9) << mesh.get_error() << std::endl;

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
  Mesh mesh(atoi(argv[1]), domain);

  /* 
  //Text Problem
  // set the function to update the mesh
  std::function<double(double, double)> f = [](double x, double y) {
    return 8*std::numbers::pi*std::numbers::pi*sin(2*std::numbers::pi*x)*cos(2*std::numbers::pi*y);
  };

  // set the exact solution
  std::function<double(double, double)> u = [](double x, double y) {
    return sin(2*std::numbers::pi*x)*cos(2*std::numbers::pi*y);
  };
  */

  
  // One possibile problem with a correct solution
  std::function<double(double, double)> f = [](double x, double y) {
    return 8*std::numbers::pi*std::numbers::pi*sin(2*std::numbers::pi*x)*sin(2*std::numbers::pi*y);
  };

  // set the exact solution
  std::function<double(double, double)> u = [](double x, double y) {
    return sin(2*std::numbers::pi*x)*sin(2*std::numbers::pi*y);
  };
  

  solution_finder(mesh, f);

  mesh.write("vtk_files/approx_sol.vtk");

  // compute the error between the mesh and the exact solution
  dist_exact_sol(mesh, u);

  return 0;
}
