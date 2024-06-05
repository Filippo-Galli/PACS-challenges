#include "Mesh.hpp"

std::optional<std::string> Mesh::parser_creation(const std::string & f){
  /**
   * @brief Function to create the parser
   * @param f is the function to parse
   * @return nullptr or an error message
  */

  p.DefineConst("pi", M_PI);
  p.DefineConst("e", M_E);

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
   * @param parser is the muparser parser
   * @return the value of the function f at x and y
  */

  parser.DefineVar("x", &x);
  parser.DefineVar("y", &y);
  return parser.Eval();
}

Mesh::Mesh(const size_t & row_number, const size_t & col_number, const Domain & domain_, const std::string & f) : mesh_data_class(row_number, col_number, domain_) {
  /**
   * @brief Constructor of the Mesh class
   * @param row_number is the number of rows of the mesh
   * @param col_number is the number of columns of the mesh
   * @param domain_ is the domain of the mesh
   * @param f is the function to be solved
  */
  f_str = f;

  auto error = parser_creation(f);
  if(error.has_value()){
    throw std::runtime_error(error.value());
  }
}

Mesh::Mesh(const std::vector<double> & _mesh, const size_t & col_number, const Domain & domain_, const std::string & f): mesh_data_class(_mesh, col_number, domain_){
  /**
   * @brief Constructor of the Mesh class
   * @param _mesh is the mesh vector
   * @param col_number is the number of columns of the mesh
   * @param domain_ is the domain of the mesh
   * @param f is the function to be solved
  */

  f_str = f;
  
  auto error = parser_creation(f);
  if(error.has_value()){
    throw std::runtime_error(error.value());
  }
}

void Mesh::update_seq(){
  /**
   * @brief Function to update the mesh using the Jacobi method sequentially
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
   * @brief Function to update the mesh using the Jacobi method - openMP parallel version
   * @param n_tasks is the number of parallel tasks
  */

  // swap the meshes, in this way the useless value are overwrite
  mesh_old.swap(mesh);

  std::pair<double, double> coords;
  #pragma omp parallel for num_threads(n_tasks) private(coords)
  for(size_t r = 1; r < n_row - 1; ++r) {
    for(size_t c = 1; c < n_col - 1; ++c) {
      // calculate + update value in the mesh
      coords = get_coordinates(r, c);
      mesh[r*n_col + c] = 0.25*(mesh_old[(r-1)*n_col + c] + mesh_old[(r+1)*n_col + c] + mesh_old[r*n_col + (c-1)] + mesh_old[r*n_col + (c+1)] + h*h*f(coords.first, coords.second, p));
    }
  }

  update_error();
}

void Mesh::update_error() { 
  /**
   * @brief Function to update the error of the current mesh with the previous one
  */

  error = 0;

  #pragma omp parallel for reduction(+:error)
  for(size_t i = 1; i < n_row - 1; ++i) 
    for(size_t j = 1; j < n_col - 1; ++j) 
      error += (mesh[i*n_col + j] - mesh_old[i*n_col + j])*(mesh[i*n_col + j] - mesh_old[i*n_col + j]);
    
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

std::optional<std::string> Mesh::add_boundary_condition(const std::string & _f){
  /**
   * @brief Function to add a boundary condition to the mesh
   * @param _f is the function of the boundary condition
  */

  mu::Parser b_parser;
  b_parser.DefineConst("pi", M_PI);
  b_parser.DefineConst("e", M_E);

  try{
    b_parser.SetExpr(_f);
  } 
  catch (mu::Parser::exception_type &e) {
    return e.GetMsg();
  }

  // add first row and last row
  for(size_t i = 0; i < n_col; ++i){
    mesh[i] = f(get_coordinates(0, i).first, get_coordinates(0, i).second, b_parser);
    mesh[(n_row - 1)*n_col + i] = f(get_coordinates(n_row - 1, i).first, get_coordinates(n_row - 1, i).second, b_parser);
  }

  // add first column and last column
  for(size_t i = 0; i < n_row; ++i){
    mesh[i*n_col] = f(get_coordinates(i, 0).first, get_coordinates(i, 0).second, b_parser);
    mesh[i*n_col + n_col - 1] = f(get_coordinates(i, n_col - 1).first, get_coordinates(i, n_col - 1).second, b_parser);
  }

  return std::nullopt;
}