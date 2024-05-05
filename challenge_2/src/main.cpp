#include "Matrix.hpp"

#include <chrono>

void prod_matrix_vector(auto &m)
{
  /**
   * @brief function to test the time taken by the product of a matrix with a vector
   * @note the function will print the result of the product and the time taken by the operation
   * @param m the matrix to test
   */
  std::vector<double> v(m.get_cols(), 1);

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<double> result = m * v;

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  
  // printing the time
  std::string placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
  std::cout << "[Matrix-Vector] " <<placeholder << duration.count() << " mus" << std::endl;

  // printing the result
  //std::cout << "Result: ";
  //for (size_t i = 0; i < result.size(); i++)
  //{
  //  std::cout << result[i] << " ";
  //}
  //std::endl(std::cout);
}

void prod_matrix_matrix(auto &m, size_t n_cols = 2)
{
  /**
   * @brief function to test the product of a matrix with another matrix
   * @note the function will print the result of the product and the time taken by the operation
   * @param m first matrix to multiply
   * @param n_cols number of columns of the second matrix
   */

  //algebra::Matrix<double, algebra::StorageOrder::ColumnMajor>m2(n_cols, m.get_cols());
  algebra::Matrix<double, algebra::StorageOrder::RowMajor> m2(m.get_cols(), n_cols);
  for(size_t i = 0; i < m.get_cols(); i++){
    if(m2.get_order() == algebra::StorageOrder::RowMajor){
      for(size_t j = 0; j < n_cols; j++){
        m2(i, j) = 1;
      }
    }
    else{
      for(size_t j = 0; j < n_cols; j++){
        m2(j, i) = 1;
      }
    } 
  }


  auto start = std::chrono::high_resolution_clock::now();

  auto m_temp = m * m2;

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  // printing the time
  std::string placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
  std::cout << "[Matrix-Matrix] "<< placeholder << duration.count() << " mus" << std::endl;

  // printing the result
  //std::cout << "Result: " << std::endl;
  //m_temp.print();
}

void norm_test(auto &m)
{
  /**
   * @brief function to test the norm of the matrix
   * @note the function will print the norm of the matrix in different norms and then compress/uncompress the matrix and print the norm again
   * @param m the matrix to test
   */

  std::string placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
  std::cout << std::endl
            << "Norm test" << std::endl
            << placeholder << std::endl;
  std::cout << "Norm - One: " << m.template norm<algebra::norm_type::One>() << std::endl;
  std::cout << "Norm - Infinity: " << m.template norm<algebra::norm_type::Infinity>() << std::endl;
  std::cout << "Norm - Frobenius: " << m.template norm<algebra::norm_type::Frobenius>() << std::endl;

  if (m.is_compressed()){
    m.uncompress();
  }
  else{
    m.compress();
  }

  placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
  std::cout << std::endl
            << placeholder << std::endl;
  std::cout << "Norm - One: " << m.template norm<algebra::norm_type::One>() << std::endl;
  std::cout << "Norm - Infinity: " << m.template norm<algebra::norm_type::Infinity>() << std::endl;
  std::cout << "Norm - Frobenius: " << m.template norm<algebra::norm_type::Frobenius>() << std::endl;
}

int main(int argc, char *argv[])
{

  if (argc == 1)
  {
    std::cout << "Usage: ./main <filename matrix 1>" << std::endl;
    return 1;
  }
  else if (argc > 2)
  {
    std::cout << "Too many arguments" << std::endl;
    return 1;
  }

  std::string filename = argv[1];

  //algebra::Matrix<double, algebra::StorageOrder::ColumnMajor> m(filename);
  algebra::Matrix<double, algebra::StorageOrder::RowMajor> m(filename);

  prod_matrix_vector(m);
  prod_matrix_matrix(m);

  m.compress();

  prod_matrix_vector(m);
  prod_matrix_matrix(m);

  //norm_test(m);


  return 0;
}
