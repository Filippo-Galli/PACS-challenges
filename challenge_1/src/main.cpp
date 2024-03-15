#include "util.hpp"
#include "config.hpp"

template<int strategy>
format decay(const parameters & p, const vector & x, const int k, const gradient_wrapper & grad, const function_wrapper & f){

  /** 
   *  @brief Function to decay the learning rate
   *  @param p: parameters
   *  @param x: vector of variables
   *  @param k: iteration number
   *  @return new learning rate
   */

  // exponential decay
  if constexpr (strategy == 1) 
    return p.alpha_0 * exp(-p.mu * k);

  // inverse decay
  else if constexpr (strategy == 2) 
    return p.alpha_0 / (1 + p.mu * k);
  
  // approximate line search with Armijo rule
  else if constexpr (strategy == 3 and mode == 0) { 
    // check condition 
    bool exit = false;

    // temporary variables
    vector temp_grad = grad(x);
    vector temp_x = subtraction_vector(x, scalar_vector(temp_grad, p.alpha_0));
    
    format temp_alpha = p.alpha_0;
    while (!exit) {

      // check condition of Amijo rule
      format temp_norm = norm2(temp_grad);
      if(f(x) - f(temp_x) >= p.sigma * temp_alpha * temp_norm * temp_norm)
        exit = true; 
      else{
        temp_alpha /= 2;
        temp_x = subtraction_vector(x, scalar_vector(temp_grad, temp_alpha));
      }
    }

    return temp_alpha;
  }
}

void strategy_chooser(format & alpha, bool & exit, bool & err, const parameters & p, 
                      const vector & x, const gradient_wrapper & grad, const function_wrapper & f, 
                      const int k = 1){
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
      alpha = decay<1>(p, x, k, grad, f);
      break;
  case 2:
      alpha = decay<2>(p, x, k, grad, f);
      break;
  default:
      if constexpr (mode == 0){
        alpha = decay<3>(p, x, k, grad, f);
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

    // Calculate the value of d (d_0)
    d = scalar_vector(temp_grad, -alpha);

    while(!exit and k < p.max_iter){
      ++k;

      // update and store value of x and old x
      x_old = x;
      x = sum_vector(x, d);

      // stopping criteria variables update
      x_diff = subtraction_vector(x, x_old);

      // check stopping criteria
      if( norm2(temp_grad) < p.residual || norm2(x_diff) < p.step_length)
        exit = true;

      else{
        // select strategy to decay alpha
        strategy_chooser(alpha, exit, err, p, x, grad, f, k);

        // update value of d
        temp_grad = grad(x); 
        temp1 = scalar_vector(temp_grad, alpha);
        temp2 = scalar_vector(d, p.nu);
        d = subtraction_vector(temp2, temp1);

      }
    }

    if(!err){
      display_result(x, temp_grad, x_diff, k, f);
    }
  }

#elif mode == 2 // Nesterov mode
  void gradient_descent(const parameters &p, const gradient_wrapper & grad, const function_wrapper & f) {
    /** 
     *  @brief Nesterov accelerated gradient descent algorithm
     *  @param p: parameters
     *  @return vector of variables representing the optimized position
     */

    // x_old = x_{k-1}
    // x = x_k
    // x_new = x_{k+1} 

    int k = 0;
    bool err = false, exit = false;
    format alpha = p.alpha_0;
    vector x_old(p.x0), temp_grad(grad(x_old));
    vector x(x_old), x_new(x.size()), x_diff(x.size()), y(x.size());

    // Calculate x_1
    x = subtraction_vector(x_old, scalar_vector(temp_grad, alpha));

    while(!exit and k < p.max_iter) {
      ++k;      

      //set stopping condition variables
      x_diff = subtraction_vector(x, x_old);
      temp_grad = grad(x);

      // Check stopping criteria
      if(norm2(temp_grad) < p.residual || norm2(x_diff) < p.step_length) {
        exit = true;
      }
      else{
        // select strategy to decay alpha
        strategy_chooser(alpha, exit, err, p, x, grad, f, k);

        // Calculate y
        y = sum_vector(x, scalar_vector(x_diff, p.nu));

        // Calculate x_new
        temp_grad = grad(y);
        x_new = subtraction_vector(y, scalar_vector(temp_grad, alpha));
        
        // Update x_old and x for the next iteration
        x_old = x;
        x = x_new;
      }
    }

    if(!err){
      display_result(x, temp_grad, x_diff, k, f);
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

    // save the old value of x
    x_old = x;

    // update x
    temp_grad = grad(x);
    x = subtraction_vector(x, scalar_vector(temp_grad, alpha));

    while(!exit and k < p.max_iter){
      ++k;

      // check variables update
      x_diff = subtraction_vector(x, x_old);
      temp_grad = grad(x);

      if(norm2(temp_grad) < p.residual || norm2(x_diff) < p.step_length)
        exit = true;
      else{
        strategy_chooser(alpha, exit, err, p, x, grad, f, k);

        // save the old value of x
        x_old = x;

        // update x
        temp_grad = grad(x);
        x = subtraction_vector(x, scalar_vector(temp_grad, alpha));
      }

    }

    if(!err){
      display_result(x, temp_grad, x_diff, k, f);
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