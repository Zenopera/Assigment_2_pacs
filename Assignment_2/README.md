Matrix Library
This is a C++ library for working with matrices, providing functionalities for matrix manipulation, compression, decompression, matrix-vector multiplication, matrix-matrix multiplication, and reading matrices from Matrix Market files.

Features
Matrix Compression: The library supports compressing matrices for efficient storage and multiplication.
Matrix Operations: Supports various matrix operations including multiplication, addition, subtraction, and norm computation.
Matrix Market Support: Can read matrices from Matrix Market files, a common format for representing sparse matrices.
Template-based Design: The library uses templates to support matrices of different data types and storage orders (Row-Major or Column-Major).
Usage
Installation
Clone the repository:
bash
Copy code
git clone https://github.com/your_username/matrix-library.git
Build the library using your preferred build system.
Example Usage
cpp
Copy code
#include "Matrix.hpp"
#include "MatrixIMPL.hpp"
#include <iostream>

int main() {
// Create a matrix
algebra::Matrix<int, algebra::StorageOrder::RowMajor> mat(3, 3);

    // Insert elements using the call operator
    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(0, 2) = 3;
    // Add more matrix operations here...

    return 0;
}



