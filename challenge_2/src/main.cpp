#include "Matrix.hpp"

int main(int argc, char *argv[]){
    /*
    algebra::Matrix<int, algebra::StorageOrder::RowMajor> m(4,6);

    // Fill the matrix
    m(0, 0) = 10;
    m(0, 1) = 20;
    m(1, 1) = 30;
    m(1, 3) = 40;
    m(2, 2) = 50;
    m(2, 3) = 60;
    m(2, 4) = 70;
    m(3, 5) = 80;

    // Print starting matrix
    std::cout << std::endl;
    m.print();
    std::cout << std::endl;

    // test compress
    m.compress();

    std::cout << std::endl;
    m.print();
    std::cout << std::endl;
    
    // test uncompress
    m.uncompress();

    std::cout << std::endl;
    m.print();
    std::cout << std::endl;

    // test resize
    m.resize(3, 4);

    std::cout << std::endl;
    m.print();
    std::cout << std::endl;

    // test operator*
    std::vector<int> v = {1, 1, 1, 1};
    std::vector<int> result = m * v;

    std::cout << "Result of the matrix-vector multiplication: ";
    for(auto & elem : result){
        std::cout << elem << ' ';
    }

    std::cout << std::endl << m(1, 1)<<std::endl;
    */
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

    m.print();

    return 0;
}
