#include "Matrix.hpp"

#include <chrono>

void time_test(auto & m){
    std::vector<double> v(m.get_cols(), 1);

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> result = m * v;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    // printing the time
    std::string placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
    std::cout << placeholder << duration.count() << " mus" << std::endl;
}

void norm_test(auto & m){
    /**
     * @brief function to test the norm of the matrix
     * @note the function will print the norm of the matrix in different norms and then compress/uncompress the matrix and print the norm again
     * @param m the matrix to test
    */

    std::string placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
    std::cout << std::endl <<"Norm test"<< std::endl << placeholder << std::endl;
    std::cout << "Norm - One: " << m.template norm<algebra::norm_type::One>() << std::endl;
    std::cout << "Norm - Infinity: " << m.template norm<algebra::norm_type::Infinity>() << std::endl;
    std::cout << "Norm - Frobenius: " << m.template norm<algebra::norm_type::Frobenius>() << std::endl;

    if(m.is_compressed()){
      m.uncompress();
    }
    else{
      m.compress();
    }

    placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
    std::cout << std::endl << placeholder << std::endl;
    std::cout << "Norm - One: " << m.template norm<algebra::norm_type::One>() << std::endl;
    std::cout << "Norm - Infinity: " << m.template norm<algebra::norm_type::Infinity>() << std::endl;
    std::cout << "Norm - Frobenius: " << m.template norm<algebra::norm_type::Frobenius>() << std::endl;
}

int main(int argc, char *argv[]){

    if(argc == 1){
        std::cout << "Usage: ./main <filename>" << std::endl;
        return 1;
    }
    else if(argc > 2){
        std::cout << "Too many arguments" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    algebra::Matrix<double, algebra::StorageOrder::ColumnMajor> m(filename);
    //algebra::Matrix<double, algebra::StorageOrder::RowMajor> m(filename);
 
    time_test(m);

    m.compress();

    time_test(m);

    norm_test(m);

    return 0;
}
