#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <array>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <fstream>
#include <sstream>

#ifdef DEBUG
#define DEBUG_MSG(msg) std::cout << msg << std::endl;
#else
#define DEBUG_MSG(msg)
#endif

namespace algebra{

    enum StorageOrder{
        RowMajor,
        ColumnMajor
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

                // check if the file is valid
                struct stat buffer;   
                if(stat(filename.c_str(), &buffer) == -1){
                    throw std::runtime_error("file not found");
                } 

                std::ifstream file(filename);
                int nonzeros = 0;

                // skip comments + read the number of rows, columns and nonzeros
                std::string line;
                bool exit = false;
                
                // TODO: try to improve it, use 1 between getline and istringstream
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
                //std::endl(std::cout);

                file.close();

            }
        
        public:

            // Friend function declaration
            template <typename U, StorageOrder Order_op>
            friend std::vector<U> operator*(Matrix<U, Order_op> & m, const std::vector<U> & v);
            
            Matrix(const size_t & idx1, const size_t & idx2){
                
                /**
                 * @brief Constructor for the Matrix class
                 * @note This constructor will ask the user for the number of rows and columns of the matrix
                 * and then will ask for the values of the matrix
                 * @param idx1 The number of index 1 [if row major, it is the number of rows, if column major, it is the number of columns]
                 * @param idx2 The number of index 2
                 * 
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

            void print() const{
                 /**
                 * @brief Print the matrix
                 * @note This function will print the matrix in the console
                */

                // TODO: adding the possibility to save to a file
                if(compressed){

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
                        std::cout << idx << " ";
                    std::cout << std::endl;

                }
                else{
                    for(size_t i = 0; i < rows; i++){
                        for(size_t j = 0; j < cols; j++){
                            if constexpr(Order == StorageOrder::RowMajor)
                                if( data.find({i, j}) != data.cend())
                                    std::cout << data.at({i, j}) << " ";
                                else  
                                    std::cout << "0 ";
                            else
                                if( data.find({j, i}) != data.cend() )
                                    std::cout << data.at({j, i}) << " ";
                                else  
                                    std::cout << "0 ";
                        }
                        std::cout << std::endl;
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

                DEBUG_MSG("Operator() non-const");
                if(check_indexes(index1, index2)){
                    if(!is_compressed()){
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
                    }
                }
                
                throw std::out_of_range("[operator()]Invalid indexes");
            
            }

            T operator()(const size_t & index1, const size_t & index2) const{
                /**
                 * @brief Get a copy of the value of the matrix
                 * @note This function will return the value of the matrix in the position (index1, index2)
                 * @param index1 The row index
                 * @param index2 The column index
                 * @return The value of the matrix in the position (index1, index2)
                */

                DEBUG_MSG("Operator() const");
                if(check_indexes(index1, index2)){
                    if(!is_compressed())
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

            void compress(){
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

                for(auto it = data.begin(); it != data.end(); it++){
                    ++compressed_data.inner_idx[it->first[0] + 1]; // add 1 to the inner index for each element in that row
                    
                    compressed_data.outer_idx.push_back(it->first[1]);
                    compressed_data.data.push_back(it->second);
                }

                // Calculate the cumulative sum so if some line is empty, the inner_idx will be correct
                int temp_sum = 0;
                for(auto & idx : compressed_data.inner_idx){
                    temp_sum += idx;
                    idx = temp_sum;
                }
                
                // Set internal state
                compressed = true;

                // map can be cleared
                data.clear();
            }

            void uncompress(){
                /**
                 * @brief Uncompress the matrix
                 * @note This function will uncompress the matrix from a CompressedMatrix struct using Compressed Sparse Row (CSR) format 
                 * or Compressed Sparse Column (CSC) format
                */

                if(!compressed)
                    return;

                data.clear();

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


            void resize(const size_t & idx1, const size_t & idx2){
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

                if(is_compressed())
                    uncompress();

                std::cout << std::endl << "[INFO] Matrix resized to " << rows << "x" << cols << std::endl;

                for(auto temp = data.begin(); temp != data.end();){
                    if(check_indexes(temp->first[0], temp->first[1]) == false){
                        temp = data.erase(temp);
                    } else {
                        ++temp;
                    }
                }

            }
                
            // Getter
            size_t get_rows() const { return rows; }
            size_t get_cols() const { return cols; }
            bool is_compressed() const { return compressed; }

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
            for(size_t i = 0; i < m.get_rows(); i++){
                for(size_t j = 0; j < m.get_cols(); j++){
                    if constexpr(Order == StorageOrder::RowMajor)
                        result[i] += m(i, j) * v[j];
                    else
                        result[j] += m(j, i) * v[i];
                }
            }

            return result;
       }
       else{
            for(size_t i = 0; i < m.compressed_data.inner_idx.size() - 1; i++){
                for(size_t j = m.compressed_data.inner_idx[i]; j < m.compressed_data.inner_idx[i + 1]; j++){
                    result[i] += m.compressed_data.data[j] * v[m.compressed_data.outer_idx[j]];
                }
            }
       }  
       return result;
    }
}

#endif
