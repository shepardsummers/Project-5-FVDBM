#include "smm.h"
#include <immintrin.h>
#include <iostream>

double* mult_mat2D_brute(double* mat1, double* mat2, int row1, int col1, int row2, int col2) {
    if (col1 != row2) {
        std::cerr << "Incompatible matrix dimensions for multiplication.\n";
        return nullptr;
    }

    double* out = new double[row1 * col2]();

    double* mat2_T = tr033(mat2, row2, col2);  // Transposed is col2 x row2
    for (int j = 0; j < row1; ++j) {
        for (int i = 0; i < col2; ++i) {
            __m256d sum = _mm256_setzero_pd();
            int k = 0;
            for(; k <= col1 - 4; k += 4) {
                __m256d vec1 = _mm256_loadu_pd(&mat1[j * col1 + k]);
                __m256d vec2 = _mm256_loadu_pd(&mat2_T[i * col1 + k]);
                sum = _mm256_fmadd_pd(vec1, vec2, sum); // fused multiply-add
                
            }
            // Horizontal sum of AVX register
            double temp[4];
            _mm256_storeu_pd(temp, sum);
            out[j*col2 + i] = temp[0] + temp[1] + temp[2] + temp[3];
            // Handle remaining elements
            for (; k < col1; ++k) {
                out[j*col2 + i] += mat1[j * col1 + k] * mat2_T[i * col1+ k];
            }

            //std::cout << "out[" << j << "][" << i << "] = " << out[j] << "\n";
        }
    }
    
    delete[] mat2_T;  // Clean up
    return out;
}

double* tr033(double* mat, int rows, int cols) {
    double* transposed = new double[rows * cols];
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            transposed[c * rows + r] = mat[r * cols + c];
        }
    }
    return transposed;
}

void emult_mat2D(double* mat_out, double* mat1, double* mat2, int size) {
    int i = 0;
    for (; i + 4 <= size; i += 4) {
        __m256d data1 = _mm256_loadu_pd(&mat1[i]); // Load 4 elements of mat1
        __m256d data2 = _mm256_loadu_pd(&mat2[i]); // Load 4 elements of mat2
        __m256d res = _mm256_mul_pd(data1, data2); // Multiply by value
        _mm256_storeu_pd(&mat_out[i], res); // Store the result
    }
    for(; i < size; ++i) {
        mat_out[i] = mat1[i]*mat2[i]; // Handle the rest of the data
    }
}

double** emult_mat2D(double* mat1, double* mat2, int size, int row) {
    double* mat_out = new double[size];
    int i = 0;
    for (; i + 4 <= size; i += 4) {
        __m256d data1 = _mm256_loadu_pd(&mat1[i]); // Load 4 elements of mat1
        __m256d data2 = _mm256_loadu_pd(&mat2[i]); // Load 4 elements of mat2
        __m256d res = _mm256_mul_pd(data1, data2); // Multiply by value
        _mm256_storeu_pd(&mat_out[i], res); // Store the result
    }
    for(; i < size; ++i) {
        mat_out[i] = mat1[i]*mat2[i]; // Handle the rest of the data
    }
    double** mat_out2D = new double*[i * row];
    for (int i = 0; i < row; ++i) {
        mat_out2D[i] = &mat_out[i * row];
    }
    return mat_out2D;
}

void scale_mat2D(double* mat, int size, double val) {
    __m256d val_vec = _mm256_set1_pd(val); // Put the value on all 4 lanes
    int i = 0;
    for (; i + 4 <= size; i += 4) {
        __m256d data = _mm256_loadu_pd(&mat[i]); // Load 4 elements
        __m256d res = _mm256_mul_pd(data, val_vec); // Multiply by value
        _mm256_storeu_pd(&mat[i], res); // Store the result
    }
    for(; i < size; ++i) {
        mat[i] *= val; // Handle the rest of the data
    }
}

void out_mat2D(double* mat, int row, int col){
    for(int j = 0; j < row; ++j) {
        for(int i = 0; i < col; ++i) {
            std::cout << mat[j*col + i] << " ";
        }
        std::cout << std::endl;
    }
}
void out_mat2D(const double* mat, int row, int col){
    for(int j = 0; j < row; ++j) {
        for(int i = 0; i < col; ++i) {
            std::cout << mat[j*col + i] << " ";
        }
        std::cout << std::endl;
    }
}

void fill_mat2D(double* p_mat, int row, int col, double val) {
    int size = row*col;
    int k = 0;
    // Put val in all 4 lanes
    __m256d vec_val = _mm256_set1_pd(val);
    // Fill in
    for(; k <= size - 4; k += 4) {
        _mm256_storeu_pd(&p_mat[k], vec_val);
    }
    // Handle remaining elements
    for (; k < size; ++k) {
        p_mat[k] = val;
    }
}