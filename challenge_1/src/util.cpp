#include <iostream>
#include <iomanip>
#include <cmath>
#include "util.hpp"

void display_parameters(const parameters& p) {
  std::cout << "\n---------- PARAMETERS ----------" << std::endl;
  std::cout << "- Max iter:" << std::setw(20) << p.max_iter << std::endl;
  std::cout << "- Step Length:" << std::setw(18) <<p.step_length << std::endl;
  std::cout << "- Residual:" << std::setw(21) << p.residual << std::endl;
  std::cout << "- Learning rate:"<< std::setw(14) << p.alpha_0 << std::endl;
  std::cout << "- Initial guess:" << std::setw(13) << "[ ";
  for(auto & i : p.x0) 
    std::cout << i << " ";

  std::cout << "]"<< std::endl;
  std::cout << "- Strategy:" << std::setw(17) << p.strategy << std::endl;
}

double norm2(const std::vector<double> & x){
  /** @brief Norm of a vector
   *  @param x: vector
   *  @return norm of the vector
   */
  double norm = 0;
  for(auto & i : x)
    norm += std::pow(i, 2);
  return std::sqrt(norm);
}

std::vector<double> scalar_vector(std::vector<double> & x, const double alpha){
  /** @brief Scalar product of a vector
   *  @param x: vector
   *  @param alpha: scalar
   *  @return scalar product
   */
  for(auto & i : x)
    i *= alpha;
  return x;
}

std::vector<double> subtraction_vector(const std::vector<double> & x1, const std::vector<double> & x2){
  /** @brief Subtraction of two vectors
   *  @param x1: vector
   *  @param x2: vector to subtract
   *  @return subtraction of the two vectors
   */
  std::vector<double> temp(x1.size());
  for(size_t i = 0; i < x1.size(); i++)
    temp[i] = x1[i] - x2[i];
  return temp;

}