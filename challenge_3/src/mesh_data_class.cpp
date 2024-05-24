#include <mesh_data_class.hpp>

#include <fstream>

mesh_data_class::mesh_data_class(const size_t & row_number, const size_t & col_number, const Domain & domain_): n_row(row_number), n_col(col_number), domain(domain_){
    /**
     * @brief Constructor of the mesh_data_class
     * @param row_number is the number of rows of the mesh
     * @param col_number is the number of columns of the mesh
    */
    
    // we use number of col since is the original n
    h = (domain.y1 - domain.y0)/(n_col-1); 

    // initialize mesh
    mesh.insert(mesh.end(), n_col*n_row, 0.0);

    // initialize mesh_old
    mesh_old = mesh;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size_th);
}

mesh_data_class::mesh_data_class(const std::vector<double> & _mesh, const size_t & col_number, const Domain & domain_): n_col(col_number), domain(domain_){
    /**
     * @brief Constructor of the mesh_data_class
     * @param _mesh is the mesh
     * @param col_number is the number of columns of the mesh
    */

    mesh = _mesh;
    
    // we use number of col since is the original n
    h = (domain.y1 - domain.y0)/(n_col-1); 

    // we use number of col since is the original n
    n_row = mesh.size()/n_col;

    // initialize mesh_old
    mesh_old = mesh;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size_th);
}

void mesh_data_class::print() const {
    /**
     * @brief Function to print the mesh
    */
    std::cout << std::endl;

    std::cout << std::setw(spacing +2) <<"| ";
    for (size_t i = 0; i < n_col; ++i) {
        std::cout << std::setw(spacing) << i << " "; // Use setw to ensure alignment
    }
    std::cout << std::endl;

    for (size_t i = 0; i < n_col + 1; ++i) {
        std::cout << std::string(spacing +1, '-');
    }
    std::cout << std::endl;

    for (size_t r = 0; r < n_row; ++r) {
      std::cout << std::setw(spacing) << r << "| "; // Use setw here as well
      for (size_t c = 0; c < n_col; ++c) {
          std::cout << std::setw(spacing) << std::setprecision (2) << mesh[r * n_col + c] << " "; // Ensure alignment for mesh values
      }
      std::cout << std::endl;
    }

    std::cout << std::endl;
}

std::optional<std::string> mesh_data_class::write(const std::string & filename) const{
  /**
   * @brief Function to write the mesh in a vtk file in Legacy format
   * @param filename is the name of the file
   * @return true if the mesh is written, false otherwise
  */
  std::ofstream file(filename);
  if (!file.is_open()) {
      return "Unable to open the file, it may not exist or you don't have the right permissions";
  }

  // Title
  file << "# vtk DataFile Version 3.0" << std::endl;
  file << "3D mesh data" << std::endl;

  // Data type
  file << "ASCII" << std::endl;

  // Geometry
  file << "DATASET STRUCTURED_GRID" << std::endl;
  file << "DIMENSIONS " << n_row << " " << n_col << " " << 1 << std::endl;

  // Points
  file << "POINTS " << n_row*n_col << " float" << std::endl;
  for(size_t i = 0; i < n_row; ++i) {
      for(size_t j = 0; j < n_col; ++j) {
        auto coords = get_coordinates(i, j);
        file << coords.first << " " << coords.second << " " << mesh[i*n_col + j] << std::endl;
      }
  }

  file.close();
  return std::nullopt;
}

std::pair<double, double> mesh_data_class::get_coordinates(const size_t & r, const size_t & c) const{

    /**
     * @brief Function to get the coordinates of the mesh
     * @param r is the row index
     * @param c is the column index
     * @return the coordinates of the mesh
    */ 

    // offset for the rank
    int offset = 0;

    if(rank == size_th - 1){
        offset = n_col - n_row;
    }
    else if (rank != 0){
        offset = rank*(n_row - 2) -1;
    }

    //std::cout << "Rank: " << rank << " r + offset: " << r + offset << std::endl;

    return std::make_pair(domain.x0 + (r + offset)*h, domain.y0 + c*h);

}

std::optional<std::string> mesh_data_class::set_mesh(const std::vector<double> & _mesh) {
    /**
     * @brief Function to set the mesh
     * @param _mesh is the mesh to be set
    */
    if(_mesh.size() != n_row*n_col){
        return "The size of the vector is not compatible with the mesh";
    }

    mesh = _mesh;
    return std::nullopt;
}

std::optional<std::string> mesh_data_class::set_mesh_old(const std::vector<double> & _mesh) {
    /**
     * @brief Function to set the mesh_old
     * @param _mesh is the mesh to be set
    */
    if(_mesh.size() != n_row*n_col){
        return "The size of the vector is not compatible with the mesh_old";
    }

    mesh_old.insert(mesh_old.end(), _mesh.begin(), _mesh.end());
    return std::nullopt;
}

