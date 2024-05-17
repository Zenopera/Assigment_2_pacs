// Created by zenop on 28/04/2024.
//
#include "Matrix.hpp"
#include "MatrixIMPL.hpp"
#include "chrono.hpp" // @note chrono is a header file that contains the Chrono class and is chrono.hpp!
#include <iostream>

// Function to perform matrix-vector multiplication and measure time
template <typename T, algebra::StorageOrder Order>
void performMultiplicationAndMeasureTime(algebra::Matrix<T, Order> &mat,
                                         const std::vector<T> &vec,
                                         const std::string &description) {
  Timings::Chrono cron;
  cron.start();

  auto res = mat * vec;

  cron.stop();

  std::cout << description << " Time: " << cron << " microseconds" << std::endl;
}

int main() {
  // Create a matrix
  algebra::Matrix<int, algebra::StorageOrder::RowMajor> mat(3, 3);

  // Insert elements using the call operator
  mat(0, 0) = 1;
  mat(0, 1) = 2;
  mat(0, 2) = 3;
  mat(1, 0) = 0;
  mat(1, 1) = 0;
  mat(1, 2) = 6;
  mat(2, 0) = 7;
  mat(2, 1) = 0;
  mat(2, 2) = 9;

  // Print the matrix
  std::cout << "Original Matrix:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::cout << mat(i, j) << " ";
    }
    std::cout << std::endl;
  }

  // Compress the matrix
  mat.compress();

  // Print the compressed state
  std::cout << "\nCompressed State:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::cout << mat(i, j) << " ";
    }
    std::cout << std::endl;
  }

  // Vector for matrix-vector multiplication
  std::vector<int> vec = {1, 2, 3};

  // Perform matrix-vector multiplication in both compressed and uncompressed
  // states
  auto result_uncompressed = mat * vec;
  mat.uncompress(); // Convert back to uncompressed state
  auto result_compressed = mat * vec;

  // Print results of matrix-vector multiplication
  std::cout << "\nMatrix-Vector Multiplication (Uncompressed):" << std::endl;
  for (const auto &elem : result_uncompressed) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;

  std::cout << "Matrix-Vector Multiplication (Compressed):" << std::endl;
  for (const auto &elem : result_compressed) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;

  // Matrix-matrix multiplication
  algebra::Matrix<int, algebra::StorageOrder::RowMajor> mat2(3, 3);
  mat2(0, 0) = 1;
  mat2(0, 1) = 2;
  mat2(0, 2) = 3;
  mat2(1, 0) = 4;
  mat2(1, 1) = 5;
  mat2(1, 2) = 6;
  mat2(2, 0) = 7;
  mat2(2, 1) = 8;
  mat2(2, 2) = 9;

  // Perform matrix-matrix multiplication in both compressed and uncompressed
  // states
  auto result_matrix_uncompressed = mat * mat2;
  mat.compress(); // Convert back to compressed state
  auto result_matrix_compressed = mat * mat2;

  // Print results of matrix-matrix multiplication
  std::cout << "\nMatrix-Matrix Multiplication (Uncompressed):" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::cout << result_matrix_uncompressed(i, j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "Matrix-Matrix Multiplication (Compressed):" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::cout << result_matrix_compressed(i, j) << " ";
    }
    std::cout << std::endl;
  }

  // Compute the Frobenius norm of the matrix
  double norm_value = mat.norm<algebra::NormType::Frobenius>();
  std::cout << "\nFrobenius Norm of the Matrix: " << norm_value << std::endl;

  // Now the matrix from matrix market:
  // Read the matrix from file
  try {
    auto m =
        algebra::Matrix<int, algebra::StorageOrder::RowMajor>::readMatrixMarket(
            "lnsp_131.mtx");
    std::vector<int> vec(m.num_Cols(),
                         1); // Generate a vector of the right dimension

    // Perform matrix-vector multiplication and measure time for different cases
    performMultiplicationAndMeasureTime(
        m, vec, "Matrix-Vector Multiplication (Uncompressed, Row-Major)");
    performMultiplicationAndMeasureTime(
        m, vec, "Matrix-Vector Multiplication (Compressed, Row-Major)");
  } catch (const std::exception &ex) {
    std::cerr << "Error: " << ex.what() << std::endl;
    return 1;
  }
  return 0;
}
