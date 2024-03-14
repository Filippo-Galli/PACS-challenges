#include "util.hpp"
#include "config.hpp"

template<int strategy>
format decay(const parameters & p, const vector & x, const int k, format alpha_0, const gradient_wrapper & grad, const function_wrapper & f){

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
  else if constexpr (strategy == 3 and mode == 0) { 
    // check condition 
    bool exit = false;

    // temporary variables
    vector temp_grad = grad(x);
    vector temp_x = subtraction_vector(x, scalar_vector(temp_grad, alpha_0));

    while (!exit) {
      // check condition of Amijo rule
      format temp_norm = norm2(temp_grad);
      if(f(x) - f(temp_x) >= p.sigma * alpha_0 * temp_norm * temp_norm)
        exit = true; 
      else{
        alpha_0 /= 2;
        temp_x = subtraction_vector(x, scalar_vector(temp_grad, alpha_0));
      }
    }

    return alpha_0;
  }
}

void strategy_chooser(format & alpha, bool & exit, bool & err,
                      const parameters & p, const vector & x, const gradient_wrapper & grad, const function_wrapper & f, const int k = 1){
  /** 
   *  @brief Function to choose the strategy
   *  @param alpha: learning rate
   *  @param exit: flag to exit the loop
   *  @param err: flag to check if there is an error
   *  @param p: parameters
   *  @param x: vector of variables  
   *  @param k: iteration number
   *  @return strategy
   */

  switch (p.strategy) {
  case 1:
      alpha = decay<1>(p, x, k, alpha, grad, f);
      break;
  case 2:
      alpha = decay<2>(p, x, k, alpha, grad, f);
      break;
  default:
      if constexpr (mode == 0){
        alpha = decay<3>(p, x, k, alpha, grad, f);
      }
      else{
        std::cerr << "\nWrong strategy" << std::endl;
        exit = true;
        err = true;
      }
      break;
    }
}

#if mode == 1 // Heavy-Ball mode
  void gradient_descent(const parameters & p, const gradient_wrapper & grad, const function_wrapper & f){
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

    vector x(p.x0), temp_grad(grad(x));
    vector x_old(x_size), x_diff(x_size), temp1(x_size), temp2(x_size), d(x_size); 

    while(!exit and k < p.max_iter){
      ++k;

      // update value of d
      temp_grad = grad(x);
      d = scalar_vector(temp_grad, -alpha);

      // update and store value of x
      x_old = x;
      x = sum_vector(x, d);
      x_diff = subtraction_vector(x, x_old);

      // set the correct gradient to check condition
      temp_grad = grad(x); 

      // check stopping criteria
      if( norm2(temp_grad) < p.residual || norm2(x_diff) < p.step_length)
        exit = true;

      else{
        // select strategy to decay alpha
        strategy_chooser(alpha, exit, err, p, x, grad, f);

        // update value of d
        temp1 = scalar_vector(temp_grad, -alpha);
        temp2 = scalar_vector(d, p.nu);
        d = sum_vector(temp2, temp1);

      }
    }

    if(!err){
      display_result(x, temp_grad, x_diff, k, f);
    }
    else{
      std::cerr << "Error in the computation" << std::endl;
    }
  }
#elif mode == 2 // Nesterov mode
  void gradient_descent(const parameters &p, const gradient_wrapper & grad, const function_wrapper & f) {
    /** 
     *  @brief Nesterov accelerated gradient descent algorithm
     *  @param p: parameters
     *  @return vector of variables representing the optimized position
     */

    int k = 0;
    bool err = false, exit = false;
    format alpha = p.alpha_0;
    vector x_old(p.x0), temp_grad(grad(x_old));
    vector x(x_old), x_diff, y;

    while(!exit and k < p.max_iter) {
      ++k;

      // Calculate x
      x = subtraction_vector(x_old, scalar_vector(temp_grad, alpha));

      // Calculate y (the lookahead position)
      x_diff = subtraction_vector(x, x_old);
      y = sum_vector(x, scalar_vector(x_diff, p.nu));

      // Calculate the gradient at the lookahead position
      temp_grad = grad(y);

      // Check stopping criteria
      if(norm2(temp_grad) < p.residual || norm2(x_diff) < p.step_length) {
        break;
      }

      // select strategy to decay alpha
      strategy_chooser(alpha, exit, err, p, x, grad, f);

      // Update x_old for the next iteration
      x_old = x;
    }

    if(!err){
      display_result(x, temp_grad, x_diff, k, f);
    }
    else{
      std::cerr << "Error in the computation" << std::endl;
    }

  }
#else // default mode 
  void gradient_descent(const parameters & p, const gradient_wrapper & grad, const function_wrapper & f){
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

    vector x_old(x_size), x_diff(x_size);
    vector x(p.x0), temp_grad(grad(x));

    while(!exit and k < p.max_iter){
      ++k;

        x_old = x;
        temp_grad = grad(x); // to set the correct gradient to check
        x = subtraction_vector(x, scalar_vector(temp_grad, alpha));
        x_diff = subtraction_vector(x, x_old);
    

      temp_grad = grad(x); // to set the correct gradient to check
        if(norm2(temp_grad) < p.residual || norm2(x_diff) < p.step_length)
          exit = true;

      

      else{
        // select strategy to decay alpha
        strategy_chooser(alpha, exit, err, p, x, grad, f);
      }
    }

    if(!err){
      display_result(x, temp_grad, x_diff, k, f);
    }
    else{
      std::cerr << "Error in the computation" << std::endl;
    }
  }
#endif

int main(){
    const parameters p;

    // display parameters
    display_parameters(p);

    // Create the function wrapper
    function_wrapper wrap_f(function);

    // Create the gradient wrapper
    #if grad_mode == 0
        gradient_wrapper wrap_grad(gradient);
    #else
        gradient_wrapper wrap_grad(function, p.h);
    #endif

    //launch the gradient descent
    gradient_descent(p, wrap_grad, wrap_f);

    return 0;
}