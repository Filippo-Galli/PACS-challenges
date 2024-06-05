#pragma once

#include "mesh_data_class.hpp"
#include <optional>

class Mesh : public mesh_data_class{
    /**
     * @brief Extension of the mesh_data_class to handle mesh update
    */

    std::optional<std::string> parser_creation(const std::string & f);
    bool check(const size_t & i, const size_t & j) const;
    double f(double x, double y, mu::Parser parser);

    public:
    Mesh(const size_t & row_number, const size_t & col_number, const Domain & domain_, const std::string & f);
    Mesh(const std::vector<double> & _mesh, const size_t & col_number, const Domain & domain_, const std::string & f); 

    // Updaters
    void update_seq();
    void update_par(const int & n_tasks = 4);
    void update_error();

    // Getters
    double get_error() const { return error; }
    std::string get_f() const { return f_str; }

    // Setters
    void set_boundary(const size_t & idx, const std::vector<double> & value, const bool & isColumn);
    std::optional<std::string> add_boundary_condition(const std::string & f);
};