//
// Created by zenop on 23/04/2024.
//
#ifndef MATRIXIMPL_H
#define MATRIXIMPL_H

#include "Matrix.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream> //@note you have forgotten this. How dould it work?
#include <sstream>


// clang-format off

namespace algebra {

    // Constructor
    template<typename T, StorageOrder Order>
    Matrix<T, Order>::Matrix(std::size_t rows, std::size_t cols)
            : numRows(rows), numCols(cols), compressed(false) {}



    // Insert element into matrix
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::insert(std::size_t row, std::size_t col, const T &value) {
        if (compressed) {
            throw std::logic_error("Cannot insert elements into a compressed matrix");
        }
        if constexpr (Order == StorageOrder::RowMajor) {
            //@note if the value exists already it will be overridden. Is this the expected behavior?
            matrixMap[{row, col}] = value;
        } else {
            matrixMap[{col, row}] = value;
        }
    }



    // Compresses the matrix representation
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::compress() {
        if (compressed) {
            //@note I think is not an error. It is just a no operation.
            std::cerr << "Already compressed" << std::endl;
            return;
        }

        // Allocate the right memory for the vectors
        outerIndexes.resize(matrixMap.size());
        values.resize(matrixMap.size());
        int k = 0;

        // Insert the values in the vectors using the right order type
        if constexpr (Order == StorageOrder::RowMajor) {
            innerIndexes.resize(numRows + 1);

            for (auto it = matrixMap.begin(); it != matrixMap.end(); it++) {
                outerIndexes[k] = it->first[1];
                values[k] = it->second;
                k++;
            }
            //@note Since map is ordered this could have been made simpler and more efficient
            // by iterating over the map and filling the vectors directly
            for (std::size_t i = 0; i < numRows + 1; i++) {
                innerIndexes[i] = std::distance(matrixMap.begin(),
                                                matrixMap.lower_bound(std::array<std::size_t, 2>{i, 0}));
            }
        } else {
            innerIndexes.resize(numCols + 1);

            for (auto it = matrixMap.begin(); it != matrixMap.end(); it++) {
                outerIndexes[k] = it->first[0];
                values[k] = it->second;
                k++;
            }

            for (std::size_t i = 0; i < numCols + 1; i++) {
                innerIndexes[i] = std::distance(matrixMap.begin(),
                                                matrixMap.lower_bound(std::array<std::size_t, 2>{0, i}));
            }
        }

        // Free the memory that contains the map
        matrixMap.clear();

        // Update the compressed status
        compressed = true;
    }



    // Uncompress matrix
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::uncompress() {
        if (!compressed) return;

        // Convert back to uncompressed format
        matrixMap.clear();
        for (std::size_t i = 0; i < numRows; ++i) {
            for (std::size_t k = innerIndexes[i]; k < innerIndexes[i + 1]; ++k) {
                if constexpr (Order == StorageOrder::RowMajor) {
                    matrixMap[{i, outerIndexes[k]}] = values[k];
                } else {
                    matrixMap[{outerIndexes[k], i}] = values[k];
                }
            }
        }

        compressed = false;
        innerIndexes.clear();
        outerIndexes.clear();
        values.clear();
    }



    // Check compression status
    template<typename T, StorageOrder Order>
    bool Matrix<T, Order>::isCompressed() const {
        return compressed;
    }



    // Resize matrix
    //@note is not just resizing, it is also clearing the matrix. 
    // This may be confusing for the user since the meaning is different from that
    // of resize for standard vectors.
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::resize(std::size_t rows, std::size_t cols) {
        numRows = rows;
        numCols = cols;
        matrixMap.clear();
        innerIndexes.clear();
        outerIndexes.clear();
        values.clear();
        compressed = false;
    }



    // Call operator overload for accessing elements (const version)
    template<typename T, StorageOrder Order>
    T Matrix<T, Order>::operator()(std::size_t row, std::size_t col) const {
        // Out of bound Check
        if constexpr (Order == StorageOrder::RowMajor) {
            if (row >= numRows || col >= numCols) {
                return T(0); // Out-of-bounds access, return 0
            }
        } else {
            if (row >= numCols || col >= numRows) {
                return T(0); // Out-of-bounds access, return 0
            }
        }

        if (compressed) {
            std::size_t index;
            if constexpr (Order == StorageOrder::RowMajor) {
                index = innerIndexes[row];
                while (index < innerIndexes[row + 1] && outerIndexes[index] < col) {
                    ++index;
                }
                if (index < innerIndexes[row + 1] && outerIndexes[index] == col) {
                    return values[index];
                }
            } else {
                index = innerIndexes[col];
             //@note if innerindexes are ordered you can use std::binary_rearch, it is more efficien
                while (index < innerIndexes[col + 1] && outerIndexes[index] < row) {
                    ++index;
                }
                if (index < innerIndexes[col + 1] && outerIndexes[index] == row) {
                    return values[index];
                }
            }
            return T(0); // Element not found, return 0
        } else {
            auto it = matrixMap.find((Order == StorageOrder::RowMajor) ? std::array{row, col} : std::array{col, row});
            if (it != matrixMap.end()) {
                return it->second;
            } else {
                return T(0); // Element not found, return 0
            }
        }
    }



    // Call operator overload for accessing elements (non-const version)
    template<typename T, StorageOrder Order>
    T &Matrix<T, Order>::operator()(std::size_t row, std::size_t col) {
        if constexpr (Order == StorageOrder::RowMajor) {
            if (row >= numRows || col >= numCols) {
                throw std::out_of_range("Row or column index out of range");
            }
        } else {
            if (row >= numCols || col >= numRows) {
                throw std::out_of_range("Row or column index out of range");
            }
        }

        if (compressed) {
            // Implement compressed matrix access (read/write)
            std::size_t index;
            if constexpr (Order == StorageOrder::RowMajor) {
                index = innerIndexes[row];
                while (index < innerIndexes[row + 1] && outerIndexes[index] < col) {
                    ++index;
                }
                if (index < innerIndexes[row + 1] && outerIndexes[index] == col) {
                    return values[index];
                } else {
                    throw std::out_of_range("Requested element is zero");
                }
            } else {
                index = innerIndexes[col];
                while (index < innerIndexes[col + 1] && outerIndexes[index] < row) {
                    ++index;
                }
                if (index < innerIndexes[col + 1] && outerIndexes[index] == row) {
                    return values[index];
                } else {
                    throw std::out_of_range("Requested element is zero");
                }
            }
        } else {
            // Uncompressed matrix access (read/write)
            return matrixMap[(Order == StorageOrder::RowMajor) ? std::array{row, col} : std::array{col, row}];
        }
    }


    // read the matrix from the file
    template<typename T, StorageOrder Order>
    Matrix<T, Order> Matrix<T, Order>::readMatrixMarket(const std::string &filename) {
        // Open the file
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        // Read the Matrix Market header
        std::string line;
        std::getline(file, line);
        if (line != "%%MatrixMarket matrix coordinate real general") {
            throw std::runtime_error("Invalid Matrix Market header");
        }

        // Read matrix size and number of non-zero elements
        std::size_t numRows, numCols, numNonZero;
        do {
            std::getline(file, line);
        } while (line[0] == '%'); // Skip comments
        std::istringstream iss(line);
        if (!(iss >> numRows >> numCols >> numNonZero)) {
            throw std::runtime_error("Invalid matrix size or non-zero count");
        }

        // Create a matrix with the read size
        Matrix<T, Order> matrix(numRows, numCols);

        // Read and insert non-zero elements
        for (std::size_t i = 0; i < numNonZero; ++i) {
            std::size_t row, col;
            T value;
            if (!(file >> row >> col >> value)) {
                throw std::runtime_error("Failed to read non-zero element");
            }
            matrix(row - 1, col - 1) = value; // Matrix Market indices are 1-based
        }

        return matrix;
    }



    template<typename T, StorageOrder Order>
    template<NormType normType>
    double Matrix<T, Order>::norm() const {
        double normValue = 0.0;

        if constexpr (Order == StorageOrder::RowMajor) {
            if constexpr (normType == NormType::One) {
                // Compute L1 norm
                // @note: In fact you can make it more efficient just by looping over the container
                // if compressied just the container holding the values, if not compressed
                // you traverse the map and extract the values. It is more efficient since you do not need
                // repeated calls to operator()(int, int). And you may also exploit standard algorithms!
                // this note applies also to the other norms.
                for (std::size_t i = 0; i < numRows; ++i) {
                    double sum = 0.0;
                    for (std::size_t j = 0; j < numCols; ++j) {
                        sum += std::abs((*this)(i, j));
                    }
                    normValue = std::max(normValue, sum);
                }
            } else if constexpr (normType == NormType::Infinity) {
                // Compute L-infinity norm
                for (std::size_t i = 0; i < numRows; ++i) {
                    double sum = 0.0;
                    for (std::size_t j = 0; j < numCols; ++j) {
                        sum += std::abs((*this)(i, j));
                    }
                    normValue = std::max(normValue, sum);
                }
            } else if constexpr (normType == NormType::Frobenius) {
                // Compute Frobenius norm
                for (std::size_t i = 0; i < numRows; ++i) {
                    for (std::size_t j = 0; j < numCols; ++j) {
                        normValue += std::pow((*this)(i, j), 2);
                    }
                }
                normValue = std::sqrt(normValue);
            }
        } else if constexpr (Order == StorageOrder::ColumnMajor) {
            if constexpr (normType == NormType::One) {
                // Compute L1 norm
                for (std::size_t j = 0; j < numCols; ++j) {
                    double sum = 0.0;
                    for (std::size_t i = 0; i < numRows; ++i) {
                        sum += std::abs((*this)(i, j));
                    }
                    normValue = std::max(normValue, sum);
                }
            } else if constexpr (normType == NormType::Infinity) {
                // Compute L-infinity norm
                for (std::size_t j = 0; j < numCols; ++j) {
                    double sum = 0.0;
                    for (std::size_t i = 0; i < numRows; ++i) {
                        sum += std::abs((*this)(i, j));
                    }
                    normValue = std::max(normValue, sum);
                }
            } else if constexpr (normType == NormType::Frobenius) {
                // Compute Frobenius norm
                for (std::size_t j = 0; j < numCols; ++j) {
                    for (std::size_t i = 0; i < numRows; ++i) {
                        normValue += std::pow((*this)(i, j), 2);
                    }
                }
                normValue = std::sqrt(normValue);
            }
        }

        return normValue;
    }


}
#endif // MATRIXIMPL_H

