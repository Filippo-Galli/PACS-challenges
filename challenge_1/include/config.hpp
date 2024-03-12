#pragma once

#include <vector>

struct parameters {
  // Strategy of decay
  // 1: exonential
  // 2: Inverse
  // 3: Approximate Line Search (Armijo rule)
  
  int strategy = 3;

  // Stopping criteria
  int max_iter = 1000;
  double step_length = 1e-6;
  double residual = 1e-6;
  
  // parameters for function decay 
  double mu = 0.2;
  double sigma = 0.1;

  // initial guess
  double alpha_0 = 0.1;
  std::vector<double> x0;

};

double f(const std::vector<double> & x);

std::vector<double> grad(const std::vector<double> & x);