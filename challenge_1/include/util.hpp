#pragma once

#include "config.hpp"

void display_parameters(const parameters& p);

double norm2(const std::vector<double> & x);

std::vector<double> scalar_vector(std::vector<double> & x, const double alpha);

std::vector<double> subtraction_vector(const std::vector<double> & x1, const std::vector<double> & x2);