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

const double* tr033(const double* mat, int rows, int cols) {
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

double* p_mult_mat2D(double* p_mat1, double* p_mat2, int row1, int col1, int page, int row2, int col2) {
    if (col1 != page) {
        std::cerr << "Incompatible matrix dimensions for multiplication.\n";
        return nullptr;
    }

    // TODO: VECTORIZE

    double* out = new double[row1*row2*col2]();
    
    for (int j = 0; j < row2; ++j) {
        for (int i = 0; i < col2; ++i) {
            for (int l = 0; l < row1; ++l) {
                for (int p = 0; p < page; ++p) {
                    //std::cout << "out[" << l << "][" << j << "][" << i << "]{" << i + (l*row2 + j)*col2 << "} = " << p_mat1[l*col1 + p] << " * " << p_mat2[i + (p*row2 + j)*col2] << std::endl;
                    out[i + (l*row2 + j)*col2]  += p_mat1[l*col1 + p]*p_mat2[i + (p*row2 + j)*col2];
                    //p + (j*col2 + i)*row2
                }
            }
        }
    }
    return out;
}

double* p_emult_mat2D(double* p_mat1, double* p_mat2, int page, int row, int col) {
    
    double* out = new double[page*row*col]();

    for (int p = 0; p < page; ++p) {
        for (int j = 0; j < row; ++j) {
            for (int i = 0; i < col; ++i) {
                out[i + (p*row + j)*col] = p_mat1[i + (p*row + j)*col]*p_mat2[i + (p*row + j)*col];
            }
        }
    }

    return out;
}

double* p_eadd_mat2D(double* p_mat1, double* p_mat2, int page, int row, int col) {
    
    double* out = new double[page*row*col]();

    for (int p = 0; p < page; ++p) {
        for (int j = 0; j < row; ++j) {
            for (int i = 0; i < col; ++i) {
                out[i + (p*row + j)*col] = p_mat1[i + (p*row + j)*col] + p_mat2[i + (p*row + j)*col];
            }
        }
    }

    return out;
}

double* p_sum_mat2D(double* p_mat, int page, int row, int col) {

    double* out = new double[row*col]();
    
    for (int j = 0; j < row; ++j) {
        for (int i = 0; i < col; ++i) {
            for (int p = 0; p < page; ++p) {
                out[i + j*col]  += p_mat[i + (p*row + j)*col];
            }
        }
    }
    return out;
}

double* p_cmult_mat2D(double* p_mat, int page, int row, int col, double val) {
    int size = page*row*col;
    double* out = new double[size]();
    for (int i = 0; i < size; ++i) {
        out[i] = p_mat[i]*val;
    }
    return out;
}

double* p_cadd_mat2D(double* p_mat, int page, int row, int col, double val) {
    int size = page*row*col;
    double* out = new double[size]();
    for (int i = 0; i < size; ++i) {
        out[i] = p_mat[i] + val;
    }
    return out;
}

double* p_emult_mat2D_a(double* p_mat1, double* p_mat2, int page, int row, int col) {
    
    double* out = new double[page*row*col]();

    for (int p = 0; p < page; ++p) {
        for (int j = 0; j < row; ++j) {
            for (int i = 0; i < col; ++i) {
                out[i + (p*row + j)*col] = p_mat1[i + (p*row + j)*col]*p_mat2[p];
            }
        }
    }

    return out;
}

double* p_eadd_mat2D_b(double* p_mat1, double* p_mat2, int page, int row, int col) {
    
    double* out = new double[page*row*col]();

    for (int p = 0; p < page; ++p) {
        for (int j = 0; j < row; ++j) {
            for (int i = 0; i < col; ++i) {
                out[i + (p*row + j)*col] = p_mat1[i + (p*row + j)*col] + p_mat2[i + j*col];
            }
        }
    }

    return out;
}

double* p_emult_mat2D_b(double* p_mat1, double* p_mat2, int page, int row, int col) {
    
    double* out = new double[page*row*col]();

    for (int p = 0; p < page; ++p) {
        for (int j = 0; j < row; ++j) {
            for (int i = 0; i < col; ++i) {
                out[i + (p*row + j)*col] = p_mat1[i + (p*row + j)*col]*p_mat2[i + j*col];
            }
        }
    }

    return out;
}