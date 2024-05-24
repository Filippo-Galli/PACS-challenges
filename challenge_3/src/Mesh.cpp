#include "Mesh.hpp"

std::optional<std::string> Mesh::parser_creation(const std::string & f){
  /**
   * @brief Function to create the parser
   * @param f is the function to parse
   * @return the parser
  */
  p.DefineConst("pi", std::numbers::pi);
  p.DefineConst("e", std::numbers::e);

  try{
    p.SetExpr(f);
  } 
  catch (mu::Parser::exception_type &e) {
    return e.GetMsg();
  }

  return std::nullopt;
}

bool Mesh::check(const size_t & i, const size_t & j) const {
    /**
     * @brief Function to check if indexes are inside the matrix 
     * @param i row index
     * @param j col index
     * @return true if indexes are inside the matrix, false otherwise  
    */
    return i < n_row && j < n_col;
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

Mesh::Mesh(const size_t & row_number, const size_t & col_number, const Domain & domain_, const std::string & f) : mesh_data_class(row_number, col_number, domain_) {
    // TODO: manage error
    parser_creation(f);

}

Mesh::Mesh(const std::vector<double> & _mesh, const size_t & col_number, const Domain & domain_, const std::string & f): mesh_data_class(_mesh, col_number, domain_){
    // TODO: manage error
    parser_creation(f);
}

void Mesh::update_seq(){
    /**
     * @brief Function to update the mesh sequentially
    */
  // swap the meshes, in this way the useless value are overwrite
  mesh_old.swap(mesh);

  for(size_t r = 1; r < n_row - 1; ++r) {
    for(size_t c = 1; c < n_col - 1; ++c) {
      // calculate + update value in the mesh
      auto coords = get_coordinates(r, c);
      mesh[r*n_col + c] = 0.25*(mesh_old[(r-1)*n_col + c] + mesh_old[(r+1)*n_col + c] + mesh_old[r*n_col + (c-1)] + mesh_old[r*n_col + (c+1)] + h*h*f(coords.first, coords.second, p));
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
  for(size_t r = 1; r < n_row - 1; ++r) {
    for(size_t c = 1; c < n_col - 1; ++c) {
      // calculate + update value in the mesh
      auto coords = get_coordinates(r, c);
      //std::cout << "Rank: " << rank << " x: " << coords.first << " y: " << coords.second << std::endl;
      mesh[r*n_col + c] = 0.25*(mesh_old[(r-1)*n_col + c] + mesh_old[(r+1)*n_col + c] + mesh_old[r*n_col + (c-1)] + mesh_old[r*n_col + (c+1)] + h*h*f(coords.first, coords.second, p));
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
    for(size_t i = 1; i < n_row - 1; ++i) {
      for(size_t j = 1; j < n_col - 1; ++j) {
        error += (mesh[i*n_col + j] - mesh_old[i*n_col + j])*(mesh[i*n_col + j] - mesh_old[i*n_col + j]);
      }
    }
    error = std::sqrt(h*error);
}


void Mesh::set_boundary(const size_t & idx, const std::vector<double> & value, const bool & isColumn){
    /**
     * @brief Function to set the boundary of the mesh
     * @param idx is the index of the boundary
     * @param value is the value of the boundary
     * @param isColumn is true if the boundary is a column, false otherwise
    */
    if(isColumn){
        for(size_t i = 0; i < n_row; ++i){
            mesh[i*n_col + idx] = value[i];
        }
    } else {
        for(size_t i = 0; i < n_col; ++i){
            mesh[idx*n_col + i] = value[i];
        }
    }
}

