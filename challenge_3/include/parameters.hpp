#pragma once

#include<functional>
#include<cmath>

struct conditions{
  double tolerance = 1e-7;
  int n_max = 1e5;

  // Distance from the exact solution for the test function 4*pi^2*cos(2*pi*x)*cos(2*pi*y) and as boundary condition 0
  // the correct solution is sin(2*pi*x)*sin(2*pi*y)
  std::function<double(double, double)> u = [](double x, double y){return sin(2*M_PI*x)*sin(2*M_PI*y);};
};