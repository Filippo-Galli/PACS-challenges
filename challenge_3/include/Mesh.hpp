#pragma once

#include<iomanip>
#include<vector>
#include<optional>
#include<iostream>
#include<functional>
#include<cmath>
#include<chrono>

struct Domain {
  double x0, x1, y0, y1;
};

class Mesh {
  std::vector<double> mesh_old;
  std::vector<double> mesh;
  size_t n;
  double h = 0;
  Domain domain;
  double error = 0;
  int spacing = 10;

  bool check(size_t i, size_t j) const;

  public:
  Mesh(const size_t & n, const Domain & domain_);

  void print(); 

  void update_seq(std::function<double(double, double)> f);

  void update_error();

  bool write(const std::string & filename) const;

  // Getters
  double get_mesh_spacing() const { return h; }
  Domain get_domain() const { return domain; }
  size_t size() const { return n; }
  std::pair<double, double> get_coordinates(size_t r, size_t c) const { return std::make_pair(domain.x0 + r*h, domain.y0 + c*h); }
  double get_error() const { return error; }
  std::optional<double> get(size_t r, size_t c) const;

  // Setters
  bool set(size_t r, size_t c, double value);
  void set_boundary(size_t idx, double value, bool isColumn);
};
