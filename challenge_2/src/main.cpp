#include "Matrix.hpp"

#include <chrono>

void time_test(auto & m){
    // Print starting matrix
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> v(m.get_cols(), 1);
    std::vector<double> result = m * v;

    // Print resulting matrix
    /*
    for(auto & i : result){
        std::cout << i << std::endl;
    }
    */

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::string placeholder = m.is_compressed() ? "Compressed: " : "Uncompressed: ";
    std::cout << placeholder << duration.count() << " ms" << std::endl;
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

    algebra::Matrix<double, algebra::StorageOrder::RowMajor> m(filename);
    
    time_test(m);

    
    std::cout << "Norm - One: " << m.norm<algebra::norm_type::One>() << std::endl;
    std::cout << "Norm - Infinity: " << m.norm<algebra::norm_type::Infinity>() << std::endl;
    std::cout << "Norm - Frobenius: " << m.norm<algebra::norm_type::Frobenius>() << std::endl;
    

    m.compress();

    time_test(m);
    
    std::cout << "Norm - One: " << m.norm<algebra::norm_type::One>() << std::endl;
    std::cout << "Norm - Infinity: " << m.norm<algebra::norm_type::Infinity>() << std::endl;
    std::cout << "Norm - Frobenius: " << m.norm<algebra::norm_type::Frobenius>() << std::endl;
    
    return 0;
}
