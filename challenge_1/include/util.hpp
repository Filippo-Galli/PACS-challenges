#pragma once

#include <iostream>
#include <iomanip>
#include <functional>

#include "config.hpp"

format norm2(const vector & x){
  /** @brief Norm of a vector
   *  @param x: vector
   *  @return norm of the vector
   */
  format norm = 0;
  for(auto & i : x)
    norm += i*i;
  return std::sqrt(norm);
}

vector scalar_vector(vector & x, const format alpha){
  /** @brief Scalar product of a vector
   *  @param x: vector
   *  @param alpha: scalar
   *  @return scalar product
   */
  for(auto & i : x)
    i *= alpha;
  return x;
}

vector subtraction_vector(const vector & x1, const vector & x2){
  /** @brief Subtraction of two vectors
   *  @param x1: vector
   *  @param x2: vector to subtract
   *  @return subtraction of the two vectors
   */
  vector temp(x1.size());
  for(size_t i = 0; i < x1.size(); i++)
    temp[i] = x1[i] - x2[i];
  return temp;

}

vector sum_vector(const vector & x1, const vector & x2){
  /** @brief Sum of two vectors
   *  @param x1: vector
   *  @param x2: vector to subtract
   *  @return element-wise sum of the two vectors
   */
  vector temp(x1.size());
  for(size_t i = 0; i < x1.size(); i++)
    temp[i] = x1[i] + x2[i];
  return temp;
}

vector grad_approx(const vector& x, std::function<format(const vector &)> f, format h = 1e-6) {
    /**
     * @brief Approximate the gradient of a function using finite differences
     * @param x: point where to compute the gradient
     * @param f: function to differentiate
     * @param h: step size
     * @return gradient of the function at x
    */

    vector grad(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        vector x_plus_h = x;
        vector x_minus_h = x;
        x_plus_h[i] += h;
        x_minus_h[i] -= h;
        grad[i] = (f(x_plus_h) - f(x_minus_h)) / (2 * h);
    }
    return grad;
}

// Function wrapper
struct function_wrapper{
    /**
     * @brief Function to be minimized wrapper
    */
    std::function<format(const vector & x)> func;

    function_wrapper(std::function<format(const vector & x)> f): func(f){}

    format operator()(const vector & x) const{
        return func(x);
    }
};

// Gradient wrapper
struct gradient_wrapper{
    /**
     * @brief Gradient of the function to be minimized wrapper
    */
    #if grad_mode == 0
        std::function<vector(const vector & x)> grad;

        gradient_wrapper(std::function<vector(const vector & x)> g): grad(g){}

        vector operator()(const vector & x) const{
          return grad(x);
        }
    #else
        std::function<format(const vector & x)> f;
        format h;

        gradient_wrapper(std::function<format(const vector & x)> g, format step = 1e-6): f(g), h(step){}

        vector operator()(const vector & x) const {
          return grad_approx(x, f, h);
        }

    #endif
};

void display_parameters(const parameters& p) {
  std::cout << "\n++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "++++++++++++ PARAMETERS ++++++++++++" << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++\n" << std::endl;
  std::cout << "- Max iter:" << std::setw(20) << p.max_iter << std::endl;
  std::cout << "- Step Length:" << std::setw(18) <<p.step_length << std::endl;
  std::cout << "- Residual:" << std::setw(21) << p.residual << std::endl;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "- Learning rate:"<< std::setw(14) << p.alpha_0 << std::endl;
  std::cout << "- mu:" << std::setw(25) << p.mu << std::endl;
  std::cout << "- nu:" << std::setw(25) << p.nu << std::endl;
  std::cout << "- sigma:" << std::setw(22) << p.sigma << std::endl;
  
  std::cout << "------------------------------------" << std::endl;
  std::cout << "- Initial guess:" << std::setw(13) << "[ ";
  for(auto & i : p.x0) 
    std::cout << i << " ";
  std::cout << "]"<< std::endl;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "- Strategy:" << std::setw(17) << p.strategy << std::endl;
  std::cout << "- Mode: ";
  if(mode == 1)
    std::cout << std::setw(29) << "Heavy Ball" << std::endl;
  else if(mode == 2)
    std::cout << std::setw(27)<<"Nesterov" << std::endl;
  else
    std::cout << std::setw(26) << "Default" << std::endl;

  std::cout << "- gradient mode: ";
  if (grad_mode == 0)
    std::cout << std::setw(17) << "Defined" << std::endl;
  else
    std::cout << std::setw(21) << "Approximate" << std::endl;
}

void display_result(const vector & x, const vector & temp_grad, const vector & x_diff, const int k, std::function<format(const vector &)> f){
  std::cout << "\n++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "++++++++++++ NERD STATS ++++++++++++" << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++\n" << std::endl;
  std::cout << "- Iteration done: " << k << std::endl;
  std::cout << "- Residual: " << norm2(temp_grad) << std::endl;
  
  std::cout << "- Step: " << norm2(x_diff) << std::endl;

  std::cout << "\n++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "++++++++++++++ RESULT ++++++++++++++" << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++\n" << std::endl;
  std::cout << "Minimum found in: " << f(x) <<std::endl;
  std::cout << "X = [ ";
  for(auto & i : x)
    std::cout << i << " ";
  std::cout << "]" << std::endl;
}