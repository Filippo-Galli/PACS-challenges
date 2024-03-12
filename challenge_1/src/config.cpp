#include "config.hpp"
#include <cmath>

double f(const std::vector<double> & x){ 
  /** @brief Function to be minimized
   *  @param x: vector of variables
   *  @return value of the function
   */

  return x[0]*x[1] + 4*std::pow(x[0], 4) + std::pow(x[1], 2) + 3*x[0];
}

std::vector<double> grad(const std::vector<double> & x){
  /** @brief Gradient of the function to be minimized
   *  @param x: vector of variables
   *  @return gradient of the function
   */
  std::vector<double> grad(2);
  grad[0] = x[1] + 16*std::pow(x[0], 3) + 3;
  grad[1] = x[0] + 2*x[1];

  return grad;
}