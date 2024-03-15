#pragma once

#include <vector>
#include <cmath>

/**
 * @brief Mode of the program
 * 0 - Default
 * 1 - Heavy Ball
 * 2 - Nesterov
*/ 
#define mode 2

/**
 * @brief Use the right gradient or approximate it
 * 0 - Use the right gradient
 * 1 - Approximate the gradient
*/
#define grad_mode 1

// Define the type of the variables, in this way is easy to change update for example to float or to Eigen::VectorXd
typedef double format;
typedef std::vector<format> vector;

struct parameters {
  /* 
  Strategy of decay:
  1 - exponential
  2 - Inverse
  3 - Approximate Line Search (Armijo rule)
  */
  int strategy = 2;

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

  // step-lentgh for the approximate gradient
  format h = 1e-6;
};


format function(const vector & x){
  /** @brief Function to be minimized
   *  @param x: vector of variables
   *  @return value of the function
   */

  return x[0]*x[1] + 4*std::pow(x[0], 4) + std::pow(x[1], 2) + 3*x[0];
}


vector gradient(const vector & x){
  /** @brief Gradient of the function to be minimized
   *  @param x: vector of variables
   *  @return gradient of the function
   */
  vector grad(x.size());
  grad[0] = x[1] + 16*std::pow(x[0], 3) + 3;
  grad[1] = x[0] + 2*x[1];

  return grad;
}