//
// Created by zenop on 23/04/2024.
//
#ifndef MATRIX_H
#define MATRIX_H

#include <map>
#include <array>
#include <vector>
#include <string>


namespace algebra {

    enum class NormType {
        One,        // L1 norm
        Infinity,   // L-infinity norm
        Frobenius   // Frobenius norm
    };

    // Enumeration for matrix storage order
    enum class StorageOrder {
        RowMajor,
        ColumnMajor
    };

    template<typename T, StorageOrder Order>
    class Matrix {
    public:
        // Constructor
        Matrix(std::size_t rows, std::size_t cols);

        // Insert element into matrix
        void insert(std::size_t row, std::size_t col, const T& value);

        // Compress matrix
        void compress();

        // Uncompress matrix
        void uncompress();

        // Check compression status
        bool isCompressed() const;

        // Resize matrix
        void resize(std::size_t rows, std::size_t cols);

        // Call operator overload for accessing elements
        T operator()(std::size_t row, std::size_t col) const;
        T& operator()(std::size_t row, std::size_t col);

        // Read matrix from Matrix Market file
        static Matrix<T, Order> readMatrixMarket(const std::string& filename);

        // Template method for computing matrix norms
        template<NormType normType>
        double norm() const;

        // Matrix-vector multiplication
        std::vector<T> operator*(const std::vector<T>& vec);

        //Matrix-Matrix multiplication
        Matrix<T, Order> operator*(const Matrix<T, Order>& other);

        //getters
        std::size_t num_Cols() const{ return numCols; };
        std::size_t num_Rows() const{ return num_Rows; };

    private:
        // Internal storage for uncompressed matrix
        std::map<std::array<std::size_t, 2>, T> matrixMap;

        // Internal storage for compressed matrix (if applicable)
        std::vector<std::size_t> innerIndexes;
        std::vector<std::size_t> outerIndexes;
        std::vector<T> values;

        // Matrix dimensions
        std::size_t numRows;
        std::size_t numCols;

        // Flag to indicate compression status
        bool compressed;
    };

    // Matrix-Matrix multiplication
    template<typename T, StorageOrder Order>
    Matrix<T, Order> Matrix<T, Order>::operator*(const Matrix<T, Order> &other) {
            static_assert(std::is_arithmetic_v<T>, "Matrix multiplication requires arithmetic types");

            if (compressed && other.compressed) {
                // Handle multiplication for compressed state
                if constexpr (Order == StorageOrder::RowMajor) {
                    if (numCols != other.numRows) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    Matrix<T, Order> result(numRows, other.numCols);
                    for (std::size_t i = 0; i < numRows; ++i) {
                        for (std::size_t j = 0; j < other.numCols; ++j) {
                            T sum = 0;
                            for (std::size_t k = innerIndexes[i]; k < innerIndexes[i + 1]; ++k) {
                                sum += values[k] * other.values[other.innerIndexes[j] + other.outerIndexes[k]];
                            }
                            result(i, j) = sum;
                        }
                    }
                    return result;
                } else if constexpr (Order == StorageOrder::ColumnMajor) {
                    if (numRows != other.numCols) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    Matrix<T, Order> result(numCols, other.numRows);
                    for (std::size_t j = 0; j < other.numRows; ++j) {
                        for (std::size_t i = 0; i < numCols; ++i) {
                            T sum = 0;
                            for (std::size_t k = innerIndexes[i]; k < innerIndexes[i + 1]; ++k) {
                                sum += values[k] * other.values[other.innerIndexes[j] + other.outerIndexes[k]];
                            }
                            result(j, i) = sum;
                        }
                    }
                    return result;
                } else {
                    throw std::logic_error("Matrix multiplication not supported for mixed storage orders");
                }
            } else {
                // Handle multiplication for uncompressed state
                if constexpr (Order == StorageOrder::RowMajor) {
                    if (numCols != other.numRows) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    Matrix<T, Order> result(numRows, other.numCols);
                    for (std::size_t i = 0; i < numRows; ++i) {
                        for (std::size_t j = 0; j < other.numCols; ++j) {
                            T sum = 0;
                            for (std::size_t k = 0; k < numCols; ++k) {
                                sum += (*this)(i, k) * other(k, j);
                            }
                            result(i, j) = sum;
                        }
                    }
                    return result;
                } else if constexpr (Order == StorageOrder::ColumnMajor) {
                    if (numRows != other.numCols) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    Matrix<T, Order> result(numCols, other.numRows);
                    for (std::size_t j = 0; j < other.numRows; ++j) {
                        for (std::size_t i = 0; i < numCols; ++i) {
                            T sum = 0;
                            for (std::size_t k = 0; k < numRows; ++k) {
                                sum += (*this)(k, i) * other(j, k);
                            }
                            result(j, i) = sum;
                        }
                    }
                    return result;
                } else {
                    throw std::logic_error("Matrix multiplication not supported for mixed storage orders");
                }
            }
        }


    // Matrix-Vector multiplication
    template<typename T, StorageOrder Order>
    std::vector<T> algebra::Matrix<T, Order>::operator*(const std::vector<T> &vec) {
            std::vector<T> result;

            // Check if the matrix is compressed
            if (compressed) {
                if constexpr (Order == StorageOrder::RowMajor) {
                    if (numCols != vec.size()) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    result.resize(numRows);
                    for (std::size_t i = 0; i < numRows; ++i) {
                        T sum = 0;
                        for (std::size_t k = innerIndexes[i]; k < innerIndexes[i + 1]; ++k) {
                            sum += values[k] * vec[outerIndexes[k]];
                        }
                        result[i] = sum;
                    }
                } else if constexpr (Order == StorageOrder::ColumnMajor) {
                    if (numRows != vec.size()) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    result.resize(numCols);
                    for (std::size_t j = 0; j < numCols; ++j) {
                        T sum = 0;
                        for (std::size_t k = innerIndexes[j]; k < innerIndexes[j + 1]; ++k) {
                            sum += values[k] * vec[outerIndexes[k]];
                        }
                        result[j] = sum;
                    }
                }
            } else {
                // Handle multiplication for uncompressed state
                if constexpr (Order == StorageOrder::RowMajor) {
                    if (numCols != vec.size()) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    result.resize(numRows);
                    for (std::size_t i = 0; i < numRows; ++i) {
                        T sum = 0;
                        for (std::size_t j = 0; j < numCols; ++j) {
                            sum += matrixMap[{i, j}] * vec[j];
                        }
                        result[i] = sum;
                    }
                } else if constexpr (Order == StorageOrder::ColumnMajor) {
                    if (numRows != vec.size()) {
                        throw std::invalid_argument("Matrix dimensions mismatch for multiplication");
                    }

                    result.resize(numCols);
                    for (std::size_t j = 0; j < numCols; ++j) {
                        T sum = 0;
                        for (std::size_t i = 0; i < numRows; ++i) {
                            sum += matrixMap[{i, j}] * vec[i];
                        }
                        result[j] = sum;
                    }
                }
            }

            return result;
        }


} // namespace algebra


#endif // MATRIX_H
