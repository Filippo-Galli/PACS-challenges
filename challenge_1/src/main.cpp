#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "config.hpp"
#include "util.hpp"

template<int strategy>
double decay(const parameters & p, 
             double alpha_0, 
             std::vector<double> & x,
             const int k = 1 
            ){

  /** 
   *  @brief Function to decay the learning rate
   *  @param p: parameters
   *  @param alpha_0: initial learning rate
   *  @param x: vector of variables
   *  @param k: iteration number
   *  @return new learning rate
   */

  if constexpr (strategy == 1) 
    return alpha_0 * exp(-p.mu * k);
  else if constexpr (strategy == 2) 
    return alpha_0 / (1 + p.mu * k);
  
  else if constexpr (strategy == 3) {  
    bool exit = false;
    std::vector<double> temp_grad = grad(x);
    std::vector<double> temp_x = subtraction_vector(x, scalar_vector(temp_grad, alpha_0));

    while (!exit) {

      if(f(x) - f(temp_x) >= p.sigma * alpha_0 * std::pow(norm2(temp_grad), 2))
        exit = true; 
      else{
        alpha_0 /= 2;
        temp_x = subtraction_vector(x, scalar_vector(temp_grad, alpha_0));
      }
    }

    return alpha_0;
  }

  else {
    static_assert(strategy >= 1 && strategy <= 3, "Unknown strategy");
    exit(1);
  }
}

std::vector<double> gradient_descent(const parameters & p){
  /** 
   *  @brief Gradient descent algorithm
   *  @param p: parameters
   *  @return vector of variables
   */

  int k = 0;
  double alpha = p.alpha_0;
  bool exit = false;
  size_t x_size = p.x0.size();
  std::vector<double> x_old(x_size), x_diff(x_size), x(p.x0);
  std::vector<double> temp_grad(x_size);

  while(k < p.max_iter && !exit){
    ++k;
    temp_grad = grad(x);

    switch (p.strategy) {
    case 1:
        alpha = decay<1>(p, alpha, x, k);
        break;
    case 2:
        alpha = decay<2>(p, alpha, x, k);
        break;
    case 3:
        alpha = decay<3>(p, alpha, x, k);
        break;
    default:
        std::cerr << "Unknown strategy" << std::endl;
      }
    
    
    x_old = x;
    x = subtraction_vector(x, scalar_vector(temp_grad, alpha));
    x_diff = subtraction_vector(x, x_old);

    if(norm2(temp_grad) < p.residual || norm2(x_diff) < p.step_length)
      exit = true;
      
  }

  std::cout << "\n---------- NERD STATS ----------" << std::endl;
  std::cout << "- Iteration done: " << k << std::endl;
  std::cout << "- Residual: " << std::sqrt(norm2(temp_grad)) << std::endl;
  std::cout << "- Step: " << std::sqrt(norm2(x_diff)) << std::endl;
  
  return x;
}

int main() {
  parameters p;
  p.x0 = {0, 0};
  display_parameters(p);

  std::vector<double> x = gradient_descent(p);
  
  std::cout << "\n---------- RESULTS ----------" << std::endl;
  std::cout << "Minimum found in: " << f(x) <<std::endl;
  std::cout << "X = [ ";
  for(auto & i : x)
    std::cout << i << " ";
  std::cout << "]" << std::endl;

  return 0;
}

