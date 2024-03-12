#pragma once

#include "config.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

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

template <typename T_input, typename T_output>
T_output norm2(const T_input & x){
  /** @brief Norm of a vector
   *  @param x: vector
   *  @return norm of the vector
   */
  format norm = 0;
  for(auto & i : x)
    norm += std::pow(i, 2);
  return std::sqrt(norm);
}

template <typename T_vec, typename T_alpha>
T_vec scalar_vector(T_vec & x, const T_alpha alpha){
  /** @brief Scalar product of a vector
   *  @param x: vector
   *  @param alpha: scalar
   *  @return scalar product
   */
  for(auto & i : x)
    i *= alpha;
  return x;
}

template <typename T>
T subtraction_vector(const T & x1, const T & x2){
  /** @brief Subtraction of two vectors
   *  @param x1: vector
   *  @param x2: vector to subtract
   *  @return subtraction of the two vectors
   */
  T temp(x1.size());
  for(size_t i = 0; i < x1.size(); i++)
    temp[i] = x1[i] - x2[i];
  return temp;

}