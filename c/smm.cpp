#include "smm.h"
#include <immintrin.h>
#include <iostream>

void emult_mat2D(double* mat_out, double* mat1, double* mat2, int size) {
    int i = 0;
    for (; i + 4 <= size; i+= 4) {
        __m256d data1 = _mm256_loadu_pd(&mat1[i]); // Load 4 elements of mat1
        __m256d data2 = _mm256_loadu_pd(&mat2[i]); // Load 4 elements of mat2
        __m256d res = _mm256_mul_pd(data1, data2); // Multiply by value
        _mm256_storeu_pd(&mat_out[i], res); // Store the result
    }
    for(; i < size; i++) {
        mat_out[i] = mat1[i]*mat2[i]; // Handle the rest of the data
    }
}

void scale_mat2D(double* mat, int size, double val) {
    __m256d val_vec = _mm256_set1_pd(val); // Put the value on all 4 lanes
    int i = 0;
    for (; i + 4 <= size; i+= 4) {
        __m256d data = _mm256_loadu_pd(&mat[i]); // Load 4 elements
        __m256d res = _mm256_mul_pd(data, val_vec); // Multiply by value
        _mm256_storeu_pd(&mat[i], res); // Store the result
    }
    for(; i < size; i++) {
        mat[i] *= val; // Handle the rest of the data
    }
}

void out_mat2D(double* mat, int col, int row){
    for(int j = 0; j < col; j++) {
        for(int i = 0; i < row; i++) {
            std::cout << mat[j*row + i] << " ";
        }
        std::cout << std::endl;
    }
}
void out_mat2D(const double* mat, int col, int row){
    for(int j = 0; j < col; j++) {
        for(int i = 0; i < row; i++) {
            std::cout << mat[j*row + i] << " ";
        }
        std::cout << std::endl;
    }
}



