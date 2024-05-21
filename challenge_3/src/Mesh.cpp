#include "Mesh.hpp"
#include <fstream>

bool Mesh::check(size_t i, size_t j) const {
    /**
     * @brief Function to check if indexes are inside the matrix 
     * @param i row index
     * @param j col index
     * @return true if indexes are inside the matrix, false otherwise  
    */
    return i < n && j < n;
}

double Mesh::f(double x, double y, mu::Parser parser){
  /**
   * @brief Function to evaluate the function f
   * @param x is the x coordinate
   * @param y is the y coordinate
   * @param p is the parser
   * @return the value of the function f
  */
  parser.DefineVar("x", &x);
  parser.DefineVar("y", &y);
  return parser.Eval();
}

Mesh::Mesh(const size_t & n, const Domain & domain_, const std::string & f) : n(n), domain(domain_) {
  // Creating possible constants
  p.DefineConst("pi", std::numbers::pi);
  p.DefineConst("e", std::numbers::e);

  try{
    p.SetExpr(f);
  } 
  catch (mu::Parser::exception_type &e) {
    std::cerr << e.GetMsg() << std::endl;
    exit(1);
  }
  
  h = (domain.x1 - domain.x0)/(n-1);
  for(size_t i = 0; i < n; ++i){
    mesh.insert(mesh.end(), n, 0.0);
    mesh_old.insert(mesh_old.end(), n, 0.0);
  }
}

std::optional<double> Mesh::get(size_t r, size_t c) const {
    /**
     * @brief Function to get the value of the mesh
     * @param r is the row index
     * @param c is the column index
     * @return the value of the mesh if the indexes are inside the matrix, std::nullopt otherwise
    */

    if(!check(r, c))
        return std::nullopt;

    return mesh.at(r*n + c);
}


bool Mesh::set(size_t r, size_t c, double value) {
    /**
     * @brief Function to set the value of the mesh
     * @param r is the row index
     * @param c is the column index
     * @param value is the value to be set
     * @return true if the value is set, false otherwise
    */
    if(check(r, c)){
        mesh[r*n + c] = value;
        return true;
    }
    return false;
}


void Mesh::set_boundary(size_t idx, double value, bool isColumn = false) {
    /**
     * @brief Function to set boundary to our problem
     * @param idx is the index of row/column of the boundary
     * @param value is the value of the boundary condition 
     * @param isColumn default is false, to use column boundary or row one 
    */

    if(isColumn and check(0, idx)){
      for(size_t i = 0; i < n; ++i){
        mesh[idx + i*n] = value;
      }
    } 
    else {
      if(check(idx, 0)){
        for(size_t i = 0; i < n; ++i){
          mesh[i + idx*n] = value;
        }
      }
    }

    // update the mesh_old
    mesh_old = mesh;
}

void Mesh::print() {
    /**
     * @brief Function to print the mesh
    */
    std::cout << std::endl;

    std::cout << std::setw(spacing +2) <<"| ";
    for (size_t i = 0; i < n; ++i) {
        std::cout << std::setw(spacing) << i << " "; // Use setw to ensure alignment
    }
    std::cout << std::endl;

    for (size_t i = 0; i < n + 1; ++i) {
        std::cout << std::string(spacing +1, '-');
    }
    std::cout << std::endl;

    for (size_t r = 0; r < n; ++r) {
        std::cout << std::setw(spacing) << r << "| "; // Use setw here as well
        for (size_t c = 0; c < n; ++c) {
            std::cout << std::setw(spacing) << std::setprecision (2) << mesh[r * n + c] << " "; // Ensure alignment for mesh values
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void Mesh::update_seq() {
  /**
   * @brief Function to update the mesh using the Jacobi method
   * @param f is the function to be solved
  */

  // swap the meshes, in this way the useless value are overwrite
  mesh_old.swap(mesh);

  for(size_t r = 1; r < n - 1; ++r) {
    for(size_t c = 1; c < n - 1; ++c) {
      // calculate + update value in the mesh
      auto coords = get_coordinates(r, c);
      mesh[r*n + c] = 0.25*(mesh_old[(r-1)*n + c] + mesh_old[(r+1)*n + c] + mesh_old[r*n + (c-1)] + mesh_old[r*n + (c+1)] + h*h*f(coords.first, coords.second, p));
    }
  }

  update_error();
}

void Mesh::update_par(const int & n_tasks) {
  /**
   * @brief Function to update the mesh using the Jacobi method
   * @param f is the function to be solved
  */

  // swap the meshes, in this way the useless value are overwrite
  mesh_old.swap(mesh);

  #pragma omp parallel for num_threads(n_tasks) schedule(static)
  for(size_t r = 1; r < n - 1; ++r) {
    for(size_t c = 1; c < n - 1; ++c) {
      // calculate + update value in the mesh
      auto coords = get_coordinates(r, c);
      mesh[r*n + c] = 0.25*(mesh_old[(r-1)*n + c] + mesh_old[(r+1)*n + c] + mesh_old[r*n + (c-1)] + mesh_old[r*n + (c+1)] + h*h*f(coords.first, coords.second, p));
    }
  }

  update_error(n_tasks);
}

void Mesh::update_error(const int & n_tasks) { 
    /**
     * @brief Function to get the error of the mesh
     * @return the error of the mesh
    */
    error = 0;
    #pragma omp parallel for num_threads(n_tasks) schedule(static) reduction(+:error)
    for(size_t i = 1; i < n - 1; ++i) {
      for(size_t j = 1; j < n - 1; ++j) {
        error += (mesh[i*n + j] - mesh_old[i*n + j])*(mesh[i*n + j] - mesh_old[i*n + j]);
      }
    }
    error = std::sqrt(h*error);
}

bool Mesh::write(const std::string & filename) const{
  /**
   * @brief Function to write the mesh in a vtk file in Legacy format
   * @param filename is the name of the file
   * @return true if the mesh is written, false otherwise
  */
  std::ofstream file(filename);
  if (!file.is_open()) {
      return false;
  }

  // Title
  file << "# vtk DataFile Version 3.0" << std::endl;
  file << "3D mesh data" << std::endl;

  // Data type
  file << "ASCII" << std::endl;

  // Geometry
  file << "DATASET STRUCTURED_GRID" << std::endl;
  file << "DIMENSIONS " << n << " " << n << " " << 1 << std::endl;

  // Points
  file << "POINTS " << n*n << " float" << std::endl;
  for(size_t i = 0; i < n; ++i) {
      for(size_t j = 0; j < n; ++j) {
        auto coords = get_coordinates(i, j);
        file << coords.first << " " << coords.second << " " << mesh[i*n + j] << std::endl;
      }
  }

  file.close();
  return true;
}