#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <array>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <numeric>
#include <execution>
#include <cmath>
#include <algorithm>

#ifdef DEBUG
#define DEBUG_MSG(msg) std::cout << msg << std::endl;
#else
#define DEBUG_MSG(msg)
#endif

// TODO: check it works for complex numbers also

namespace algebra{

    enum StorageOrder{
        RowMajor,
        ColumnMajor
    };

    enum norm_type{
        One,
        Infinity,
        Frobenius,
    };

    template <typename T>
    struct CompressedMatrix{
        std::vector<size_t> inner_idx;
        std::vector<size_t> outer_idx;
        std::vector<T> data;
    };

    template <typename T, StorageOrder Order>
    class Matrix {
        private: 
            size_t rows = 0, cols = 0;
            std::map<std::array<std::size_t, 2>, T> data;

            bool compressed = false;

            CompressedMatrix<T> compressed_data;

            bool check_indexes(const size_t & index1, const size_t & index2) const{
                /**
                 * @brief Check if the indexes are valid
                 * @note This function will check if the indexes are valid
                 * @param row The row index
                 * @param col The column index
                 * @return True if the indexes are valid, false otherwise
                */
                if constexpr(Order == StorageOrder::RowMajor)
                    return index1 < rows && index2 < cols;
                else
                    return index1 < cols && index2 < rows;
            }

            void read_matrix_MM(const std::string & filename){
                /**
                 * @brief Read a matrix in Matrix Market format
                 * @note This function will read a matrix in Matrix Market format
                 * @param filename The name of the file
                */

                // check if the filename is valid
                struct stat buffer;
                if(stat(filename.c_str(), &buffer) == -1){
                    throw std::runtime_error("file not found");
                }

                std::ifstream file(filename);
                int nonzeros = 0;

                // skip comments + read the number of rows, columns and nonzeros
                std::string line;
                bool exit = false;

                while (!exit)
                {
                    getline(file, line);
                    if (line[0] != '%'){
                        std::istringstream useful_data(line);
                        useful_data >> rows >> cols >> nonzeros;
                        exit = true;
                    }
                }
                
                // read the matrix
                size_t row, col;
                T value;
                for (int i = 0; i < nonzeros; ++i) {
                    file >> row >> col >> value;
                    if constexpr(Order == StorageOrder::RowMajor)
                        data[{row - 1, col - 1}] = value;
                    else
                        data[{col - 1, row - 1}] = value;
                }

                file.close();

            }
        
        public:

            // Friend function declaration
            template <typename U, StorageOrder Order_op>
            friend std::vector<U> operator*(Matrix<U, Order_op> & m, const std::vector<U> & v);
            
            Matrix(const size_t & idx1, const size_t & idx2) noexcept{
                /**
                 * @brief Constructor for the Matrix class
                 * @note This constructor will ask the user for the number of rows and columns of the matrix
                 * and then will ask for the values of the matrix
                 * @param idx1 The number of index 1 [if row major, it is the number of rows, if column major, it is the number of columns]
                 * @param idx2 The number of index 2
                */

                if constexpr(Order == StorageOrder::RowMajor){
                    this->rows = idx1;
                    this->cols = idx2;
                }
                else{
                    this->cols = idx1;
                    this->rows = idx2;
                }
            
            }

            Matrix(const std::string & filename){
                /**
                 * @brief Constructor for the Matrix class
                 * @note This constructor will read a matrix in Matrix Market format
                 * @param filename The name of the file
                */

                read_matrix_MM(filename);
            }

            void print() const noexcept{
                 /**
                 * @brief Print the matrix
                 * @note This function will print the matrix in the console
                */

                if(compressed){
                    std::string placeholder = Order == StorageOrder::RowMajor ? "Row Major" : "Column Major";
                    std::cout << placeholder << std::endl;

                    std::cout << "Inner indexes: " << std::endl;
                    for(const auto & idx : compressed_data.inner_idx)
                        std::cout << idx << " ";
                    std::cout << std::endl;

                    std::cout << "Outer indexes: " << std::endl;
                    for(const auto & idx : compressed_data.outer_idx)
                        std::cout << idx << " ";
                    std::cout << std::endl;

                    std::cout << "Data: " << std::endl;
                    for(const auto & idx : compressed_data.data)
                        std::cout << idx << " " << std::endl;
                    std::cout << std::endl;
                }
                else{
                    // print guidelines on how I print the matrix
                    std::string placeholder1 = Order == StorageOrder::RowMajor ? "Row" : "Column";
                    std::string placeholder2 = Order == StorageOrder::RowMajor ? "Column" : "Row";
                    std::cout << placeholder1 << " " << placeholder2 << " " << " Data " << std::endl;
                    std::cout << "----------------------" << std::endl;
                    for(auto && pair: data){
                        std::cout << pair.first[0] << " " << pair.first[1] << " " << pair.second << std::endl;
                    }
                }

            }

            T & operator()(const size_t & index1, const size_t & index2){
                /**
                 * @brief Get the reference to the value of the matrix
                 * @note This function will return the value of the matrix in the position (index1, index2)
                 * @param index1 The row index
                 * @param index2 The column index
                 * @return The value of the matrix in the position (index1, index2)
                */

                //DEBUG_MSG("Operator() non-const");
                if(check_indexes(index1, index2)){
                    if(!is_compressed()){
                        // if the element is not in the map, add it
                        if(data.find({index1, index2}) == data.end())
                            data[{index1, index2}] = 0;
                
                        return data[{index1, index2}];
                    }
                    else{
                        size_t idx = compressed_data.inner_idx[index1];

                        while(idx < compressed_data.inner_idx[index1 + 1]){
                            if(compressed_data.outer_idx[idx] == index2){
                                return compressed_data.data[idx];
                            }

                            ++idx;
                        }

                        throw std::out_of_range("[operator()] Attempt to add value while the matrix is compressed");
                    }
                }
                
                throw std::out_of_range("[operator()] Invalid indexes");
            }

            T operator()(const size_t & index1, const size_t & index2) const{
                /**
                 * @brief Get a copy of the value of the matrix
                 * @note This function will return the value of the matrix in the position (index1, index2)
                 * @param index1 The row index
                 * @param index2 The column index
                 * @return The value of the matrix in the position (index1, index2)
                */

                //DEBUG_MSG("Operator() const");
                if(check_indexes(index1, index2)){
                    if(!is_compressed())
                        // if the element is not in the map, return 0
                        if(data.find({index1, index2}) == data.cend())
                            return 0;
                        else
                            return data.at({index1, index2});
                    else{
                        size_t idx = compressed_data.inner_idx[index1];
                        while(idx < compressed_data.inner_idx[index1 + 1]){
                            if(compressed_data.outer_idx[idx] == index2)
                                return compressed_data.data[idx];
                            ++idx;
                        }
                        return 0;
                    }
                }
                else{
                    throw std::out_of_range("[operator() const] Invalid indexes");
                }
            }

            void compress() noexcept{
                /**
                 * @brief Compress the matrix
                 * @note This function will compress the matrix in a CompressedMatrix struct using Compressed Sparse Row (CSR) format 
                 * or Compressed Sparse Column (CSC) format
                */

                if(compressed)
                    return;

                // allocate the memory for the compressed matrix
                if constexpr(Order == StorageOrder::RowMajor)
                    compressed_data.inner_idx.resize(rows +1, 0);
                else
                    compressed_data.inner_idx.resize(cols +1, 0);

                compressed_data.outer_idx.reserve(data.size());
                compressed_data.data.reserve(data.size());

                // temporary object to store the number of elements in each row
                auto temp(compressed_data.inner_idx);

                // fill the compressed matrix struct
                for(auto it = data.begin(); it != data.end(); it++){
                    ++temp[it->first[0] + 1]; // add 1 to the inner index for each element in that row
                    compressed_data.outer_idx.emplace_back(it->first[1]);
                    compressed_data.data.emplace_back(it->second);
                }

                // Calculate the cumulative sum so if some line is empty, the inner_idx will be correct
                std::partial_sum(temp.begin(), temp.end(), compressed_data.inner_idx.begin());
                
                // Set internal state
                compressed = true;

                // map can be cleared
                data.clear();
            }

            void uncompress() noexcept{
                /**
                 * @brief Uncompress the matrix
                 * @note This function will uncompress the matrix from a CompressedMatrix struct using Compressed Sparse Row (CSR) format 
                 * or Compressed Sparse Column (CSC) format
                */

                if(!compressed)
                    return;

                size_t idx = 0;
                for(size_t i = 0; i < compressed_data.inner_idx.size() -1; i++){
                    for(size_t j = compressed_data.inner_idx[i]; j < compressed_data.inner_idx[i + 1]; j++){
                        data[{i, compressed_data.outer_idx[j]}] = compressed_data.data[idx];
                        ++idx;
                    }
                }

                // Set internal state
                compressed = false;

                // compressed data can be cleared
                compressed_data.inner_idx.clear();
                compressed_data.outer_idx.clear();
                compressed_data.data.clear();
            }

            void resize(const size_t & idx1, const size_t & idx2) noexcept{
                /**
                 * @brief Resize the matrix
                 * @note This function will resize the matrix
                 * @param idx1 The number of index 1 [if row major, it is the number of rows, if column major, it is the number of columns]
                 * @param idx2 The number of index 2
                */
                if constexpr(Order == StorageOrder::RowMajor){
                    rows = idx1;
                    cols = idx2;
                }
                else{
                    cols = idx1;
                    rows = idx2;
                }
                
                // Check if the matrix is compressed and uncompress it
                if(is_compressed())
                    uncompress();

                std::cout << std::endl << "[INFO] Matrix resized to " << rows << "x" << cols << std::endl;

                for(auto temp = data.begin(); temp != data.end();){
                    // check if the indexes are valid with new dimensions
                    if(check_indexes(temp->first[0], temp->first[1]) == false){
                        temp = data.erase(temp);
                    } else {
                        ++temp;
                    }
                }

            }
            
            template <norm_type Norm>
            double norm() const{
                /**
                 * @brief Calculate the norm of the matrix
                 * @note This function will calculate the norm of the matrix
                 * @return The norm of the matrix
                */
                if constexpr(Norm == norm_type::One){
                    double max = 0;
                    if(!is_compressed()){
                        for(size_t c = 0; c < cols; c++){
                            T sum = 0; // temporary variable to store the sum of the elements of each row
                            for(size_t r = 0; r < rows; r++){
                                if constexpr(Order == StorageOrder::RowMajor)
                                    sum += std::abs((*this)(r, c));
                                else
                                    sum += std::abs((*this)(c, r));
                            }

                            if(sum > max)
                                max = sum;
                        }
                    }
                    else{
                        // TODO: improve version for RowMajor 
                        if constexpr (Order == StorageOrder::RowMajor){
                            for(size_t i = 0; i < cols; i++){
                                T sum = 0;
                                //search inside the outer index array if we reach the correct column
                                for(size_t j = 0; j < compressed_data.outer_idx.size(); j++){
                                    if(compressed_data.outer_idx[j] == i)
                                        sum += std::abs(compressed_data.data[j]);
                                }
                                if(sum > max)
                                    max = sum;
                            }
                        }
                        else{
                            for(size_t i = 0; i < cols; i++){
                                T sum = 0;
                                for(size_t j = compressed_data.inner_idx[i]; j < compressed_data.inner_idx[i + 1]; j++){
                                    sum += std::abs(compressed_data.data[j]);
                                }
                                
                                if(sum > max)
                                    max = sum;
                            }
                        }
                    }

                    return max;
                }
                else if constexpr(Norm == norm_type::Infinity){
                    T max = 0;
                    if(!is_compressed()){
                        for(size_t r = 0; r < rows; r++){
                            T sum = 0; // temporary variable to store the sum of the elements of each row
                            for(size_t c = 0; c < cols; c++){
                                if constexpr(Order == StorageOrder::RowMajor)
                                    sum += std::abs((*this)(r, c));
                                else
                                    sum += std::abs((*this)(c, r));
                            }

                            if(sum > max)
                                max = sum;
                        }
                    }
                    else{
                        if constexpr (Order == StorageOrder::RowMajor){
                            for(size_t i = 0; i < rows; i++){
                                T sum = 0;
                                for(size_t j = compressed_data.inner_idx[i]; j < compressed_data.inner_idx[i + 1]; j++){
                                    sum += std::abs(compressed_data.data[j]);
                                }
                                if(sum > max)
                                    max = sum;
                            }
                        }
                        // TODO: improve version for ColumnMajor
                        else{
                            for(size_t i = 0; i < rows; i++){
                                T sum = 0;
                                // search inside the outer index array if we reach the correct row
                                for(size_t j = 0; j < compressed_data.outer_idx.size(); j++){
                                  if(compressed_data.outer_idx[j] == i)
                                    sum += std::abs(compressed_data.data[j]);
                                }
                                if(sum > max)
                                    max = sum;
                            }
                        }
                    }

                    return max;
                }
                else if constexpr(Norm == norm_type::Frobenius){
                    double sum = 0;
                    if(!is_compressed()){
                        for(size_t r = 0; r < rows; r++){
                            for(size_t c = 0; c < cols; c++){
                                if constexpr(Order == StorageOrder::RowMajor)
                                    sum += (*this)(r, c) * (*this)(r, c);
                                else
                                    sum += (*this)(c, r) * (*this)(c, r);
                            }
                        }
                    }
                    else{
                        // temporary object to store the squared absolute values
                        auto temp(compressed_data.data);

                        // Transform the vector to hold the squared absolute values
                        std::transform(temp.begin(), temp.end(), temp.begin(), [](T val) { return std::abs(val) * std::abs(val); });

                        // Sum the elements of the vector
                        sum = std::accumulate(temp.begin(), temp.end(), 0.0);
                    }
                    return std::sqrt(sum);
                }

                throw std::invalid_argument("[norm] Invalid norm type");
            } 

            // Getter
            size_t get_rows() const noexcept { return rows; }
            size_t get_cols() const noexcept{ return cols; }
            bool is_compressed() const noexcept { return compressed; }

    };

    template <typename T, StorageOrder Order>
    std::vector<T> operator*(Matrix<T, Order> & m, const std::vector<T> & v){
        /**
         * @brief Multiply a matrix by a vector
         * @note This function will multiply a matrix by a vector
         * @param m The matrix
         * @param v The vector
         * @return The result of the multiplication
        */

       std::vector<T> result;
        
      // Error handling of wrong dimensions
       if constexpr(Order == StorageOrder::RowMajor){
            if(m.get_cols() != v.size()){
                throw std::invalid_argument("[Operator*] The number of columns of the matrix must be equal to the size of the vector");
            }

            result.resize(m.get_rows(), 0);
       }
       else{
            if(m.get_rows() != v.size()){
                throw std::invalid_argument("[Operator*] The number of columns of the matrix must be equal to the size of the vector");
            }

            result.resize(m.get_cols(), 0);
       }

       if(!m.is_compressed()){
            size_t temp_rows = m.get_rows(), temp_cols = m.get_cols();
            for(size_t i = 0; i < temp_rows; i++){
                for(size_t j = 0; j < temp_cols; j++){
                    if constexpr(Order == StorageOrder::RowMajor)
                        result[i] += m(i, j) * v[j];
                    else
                        result[j] += m(j, i) * v[i];
                }
            }

            return result;
       }
       else{
            std::vector<T> temp;
            for(size_t i = 0; i < m.compressed_data.inner_idx.size() - 1; i++){
                temp.resize(m.compressed_data.inner_idx[i + 1] - m.compressed_data.inner_idx[i]);
                
                // multiply the values of the compressed matrix with the values of the vector
                
                for(size_t j = m.compressed_data.inner_idx[i]; j < m.compressed_data.inner_idx[i + 1]; j++){
                    temp[i] = m.compressed_data.data[j] * v[m.compressed_data.outer_idx[j]];
                }
                
                //std::transform(std::execution::par, m.compressed_data.data.begin() + m.compressed_data.inner_idx[i], m.compressed_data.data.begin() + m.compressed_data.inner_idx[i + 1], temp.begin(), [&v, &m, i](T val) { return val * v[m.compressed_data.outer_idx[i]]; });

                // sum the values of the temp vector
                result[i] = std::reduce(std::execution::par, temp.begin(), temp.end(), 0);
            }
       }  
       return result;
    }
}

#endif
