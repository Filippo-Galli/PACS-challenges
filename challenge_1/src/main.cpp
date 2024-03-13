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

  // exponential decay
  if constexpr (strategy == 1) 
    return alpha_0 * exp(-p.mu * k);

  // inverse decay
  else if constexpr (strategy == 2) 
    return alpha_0 / (1 + p.mu * k);
  
  // approximate line search with Armijo rule
  else if constexpr (strategy == 3) { 
    // check condition 
    bool exit = false;

    // temporary variables
    T temp_grad = grad<T>(x);
    T temp_x = subtraction_vector<T>(x, scalar_vector<T, format>(temp_grad, alpha_0));

    while (!exit) {
      // check condition
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


#if mode == 1 // Heavy-Ball mode
  template<typename T>
  void gradient_descent(const parameters & p){
    /** 
     *  @brief Gradient descent algorithm with heavy-ball method
     *  @param p: parameters
     *  @return vector of variables
     */

    int k = 0;
    bool err = false;
    format alpha = p.alpha_0;
    bool exit = false;
    size_t x_size = p.x0.size();

    T x(p.x0), temp_grad(grad<T>(x));
    T x_old(x_size), x_diff(x_size), temp1(x_size), temp2(x_size), d(x_size); 
    while(!exit and k < p.max_iter){
      ++k;

      // update value of d
      temp_grad = grad<T>(x);
      d = scalar_vector<T, format>(temp_grad, -alpha);

      // update and store value of x
      x_old = x;
      x = sum_vector<T>(x, d);
      x_diff = subtraction_vector<T>(x, x_old);

      // set the correct gradient to check condition
      temp_grad = grad<T>(x); 

      // check stopping criteria
      if( norm2<T, format>(temp_grad) < p.residual || norm2<T, format>(x_diff) < p.step_length)
        exit = true;

      else{
        // select strategy to decay alpha
        switch (p.strategy) {
        case 1:
            alpha = decay<1, T>(p, alpha, x, k);
            break;
        case 2:
            alpha = decay<2, T>(p, alpha, x, k);
            break;
        default:
            std::cerr << "\nWrong strategy" << std::endl;
            exit = true;
            err = true;
          }

        // update value of d
        temp1 = scalar_vector<T, format>(temp_grad, -alpha);
        temp2 = scalar_vector<T, format>(d, p.nu);
        d = sum_vector<T>(temp2, temp1);

      }
    }

    if(!err){
      display_result<T, format>(x, temp_grad, x_diff, k);
    }
  }
#elif mode == 2 // Nesterov mode

  template<typename T>
  void gradient_descent(const parameters &p) {
    /** 
     *  @brief Nesterov accelerated gradient descent algorithm
     *  @param p: parameters
     *  @return vector of variables representing the optimized position
     */

    int k = 0;
    bool err = false, exit = false;
    format alpha = p.alpha_0;
    T x_old(p.x0), temp_grad(grad<T>(x_old));
    T x(x_old), x_diff, y;

    while(!exit and k < p.max_iter) {
      ++k;

      // Calculate x
      x = subtraction_vector<T>(x_old, scalar_vector<T, format>(temp_grad, alpha));

      // Calculate y (the lookahead position)
      x_diff = subtraction_vector<T>(x, x_old);
      y = sum_vector<T>(x, scalar_vector<T, format>(x_diff, p.nu));

      // Calculate the gradient at the lookahead position
      temp_grad = grad<T>(y);

      // Check stopping criteria
      if(norm2<T, format>(temp_grad) < p.residual || norm2<T, format>(x_diff) < p.step_length) {
        break;
      }

      // Update the learning rate based on the chosen strategy
      switch (p.strategy) {
        case 1:
          alpha = decay<1, T>(p, alpha, x, k);
          break;
        case 2:
          alpha = decay<2, T>(p, alpha, x, k);
          break;
        default:
          std::cerr << "\nWrong strategy" << std::endl;
          exit = true;
          err = true;
      }

      // Update x_old for the next iteration
      x_old = x;
    }

    if(!err){
      display_result<T, format>(x, temp_grad, x_diff, k);
    }
    
  }

#else

  template<typename T>
  void gradient_descent(const parameters & p){
    /** 
     *  @brief Gradient descent algorithm
     *  @param p: parameters
     *  @return vector of variables
     */

    int k = 0;
    bool err = false;
    format alpha = p.alpha_0;
    bool exit = false;
    size_t x_size = p.x0.size();

    T x_old(x_size), x_diff(x_size);
    T x(p.x0), temp_grad(grad<T>(x));

    while(!exit and k < p.max_iter){
      ++k;

        x_old = x;
        temp_grad = grad<T>(x); // to set the correct gradient to check
        x = subtraction_vector<T>(x, scalar_vector<T, format>(temp_grad, alpha));
        x_diff = subtraction_vector<T>(x, x_old);
    

      temp_grad = grad<T>(x); // to set the correct gradient to check
        if(norm2<T, format>(temp_grad) < p.residual || norm2<T, format>(x_diff) < p.step_length)
          exit = true;

      

      else{
        switch (p.strategy) {
        case 1:
            alpha = decay<1, T>(p, alpha, x, k);
            break;
        case 2:
            alpha = decay<2, T>(p, alpha, x, k);
            break;
        case 3:
            if constexpr (mode == 1){
              std::cerr << "\nWrong strategy with heavy-ball" << std::endl;
              exit = true;
              err = true;
            }  
            if(!err)
              alpha = decay<3, T>(p, alpha, x, k);
            break;
        default:
            std::cerr << "\nUnknown strategy" << std::endl;
            exit = true;
            err = true;
          }
      }
    }

    if(!err){
      display_result<T, format>(x, temp_grad, x_diff, k);
    }
  }

#endif

int main() {
  parameters p;
  display_parameters(p);

  gradient_descent<vector>(p);

  return 0;
}

