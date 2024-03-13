#pragma once

#include <vector>
#include <cmath>

/**
 * @brief Mode of the program
 * 0 - Default
 * 1 - Heavy Ball
 * 2 - Nesterov
*/ 
#define mode 0

// Define the type of the variables
typedef double format;
typedef std::vector<format> vector;

struct parameters {
  /* 
  Strategy of decay:
  1 - exponential
  2 - Inverse
  3 - Approximate Line Search (Armijo rule)
  */
  int strategy = 1;

  // Stopping criteria
  int max_iter = 1000;
  format step_length = 1e-6;
  format residual = 1e-6;
  
  // parameters for function decay 
  format mu = 0.2;
  format sigma = 0.1;

  // parameters for heavy_ball method
  format nu = 0.9;

  // initial guess
  format alpha_0 = 0.1;
  vector x0{0, 0};

};

/**
  *  Function to be minimized
  *  f(x) = x1*x2 + 4*x1^4 + x2^2 + 3*x1
  *  Gradient of the function
  *  grad(f) = [x2 + 16*x1^3 + 3, x1 + 2*x2]
*/
template<typename T_input, typename T_output>
T_output f(const T_input & x){
  /** @brief Function to be minimized
   *  @param x: vector of variables
   *  @return value of the function
   */

  return x[0]*x[1] + 4*std::pow(x[0], 4) + std::pow(x[1], 2) + 3*x[0];
}

template <typename T>
T grad(const T & x){
  /** @brief Gradient of the function to be minimized
   *  @param x: vector of variables
   *  @return gradient of the function
   */
  T grad(2);
  grad[0] = x[1] + 16*std::pow(x[0], 3) + 3;
  grad[1] = x[0] + 2*x[1];

  return grad;
}