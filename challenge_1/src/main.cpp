#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "config.hpp"
#include "util.hpp"

template<int strategy, typename T>
format decay(const parameters & p, 
             format alpha_0, 
             T & x,
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
    T temp_grad = grad<T>(x);
    T temp_x = subtraction_vector<T>(x, scalar_vector<T, format>(temp_grad, alpha_0));

    while (!exit) {

      if(f<T, format>(x) - f<T, format>(temp_x) >= p.sigma * alpha_0 * std::pow(norm2<T, format>(temp_grad), 2))
        exit = true; 
      else{
        alpha_0 /= 2;
        temp_x = subtraction_vector<T>(x, scalar_vector<T, format>(temp_grad, alpha_0));
      }
    }

    return alpha_0;
  }

  else {
    static_assert(strategy >= 1 && strategy <= 3, "Unknown strategy");
    exit(1);
  }
}

template<typename T>
T gradient_descent(const parameters & p){
  /** 
   *  @brief Gradient descent algorithm
   *  @param p: parameters
   *  @return vector of variables
   */

  int k = 0;
  format alpha = p.alpha_0;
  bool exit = false;
  size_t x_size = p.x0.size();
  T x_old(x_size), x_diff(x_size), x(p.x0);
  T temp_grad(x_size);

  while(k < p.max_iter && !exit){
    ++k;
    temp_grad = grad<T>(x);

    switch (p.strategy) {
    case 1:
        alpha = decay<1, T>(p, alpha, x, k);
        break;
    case 2:
        alpha = decay<2, T>(p, alpha, x, k);
        break;
    case 3:
        alpha = decay<3, T>(p, alpha, x, k);
        break;
    default:
        std::cerr << "Unknown strategy" << std::endl;
      }
    
    
    x_old = x;
    x = subtraction_vector<T>(x, scalar_vector<T, format>(temp_grad, alpha));
    x_diff = subtraction_vector<T>(x, x_old);

    if(norm2<T, format>(temp_grad) < p.residual || norm2<T, format>(x_diff) < p.step_length)
      exit = true;
      
  }

  std::cout << "\n---------- NERD STATS ----------" << std::endl;
  std::cout << "- Iteration done: " << k << std::endl;
  std::cout << "- Residual: " << std::sqrt(norm2<T, format>(temp_grad)) << std::endl;
  std::cout << "- Step: " << std::sqrt(norm2<T, format>(x_diff)) << std::endl;
  
  return x;
}

int main() {
  parameters p;
  display_parameters(p);

  vector x = gradient_descent<vector>(p);
  
  std::cout << "\n---------- RESULTS ----------" << std::endl;
  std::cout << "Minimum found in: " << f<vector, format>(x) <<std::endl;
  std::cout << "X = [ ";
  for(auto & i : x)
    std::cout << i << " ";
  std::cout << "]" << std::endl;

  return 0;
}

