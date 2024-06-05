#pragma once

#include<iomanip>
#include<vector>
#include<optional>
#include<iostream>
#include<functional>
#include<cmath>
#include<chrono>
#include<muParser.h>
#include<omp.h>
#include<mpi.h>

struct Domain {
    double x0, x1, y0, y1;
};

class mesh_data_class {
    /**
     * @brief Base class to handle the mesh data
    */

    protected:

    int spacing = 10;
    std::vector<double> mesh;
    std::vector<double> mesh_old;
    size_t n_row, n_col;
    double h;
    Domain domain;

    // MPI variables
    int rank = 0;
    int size_th = 1;

    // offset to corrected coordinates
    mutable int offset = 0;

    // muParser variables
    mu::Parser p;
    std::string f_str;

    // error
    double error = 0;
    
    public:
    mesh_data_class(const size_t & row_number, const size_t & col_number, const Domain & domain_);
    mesh_data_class(const std::vector<double> & _mesh, const size_t & col_number, const Domain & domain_);

    void print() const;

    std::optional<std::string> write(const std::string & filename) const;

    // Getters
    std::vector<double> & get_mesh() { return mesh; }
    std::vector<double> get_mesh_old() const { return mesh_old; }
    std::pair<size_t, size_t> get_size() const { return std::make_pair(n_row, n_col); }
    double get_h() const { return h; }
    Domain get_domain() const { return domain; }
    std::pair<double, double> get_coordinates(const size_t & r, const size_t & c) const;
    double get_value(const size_t & r, const size_t & c) const {return mesh[r*n_col + c];};

    // Setters
    std::optional<std::string> set_mesh(const std::vector<double> & _mesh);
    std::optional<std::string> set_mesh_old(const std::vector<double> & _mesh);
    void set_offset(const int & _offset) { offset = _offset; }
};