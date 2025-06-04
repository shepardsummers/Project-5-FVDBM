#ifndef SMM_H
#define SMM_H

/*
* Multiplies two 2D matricies. Columns of first must match rows of second.
* @param mat_out Pointer to first element of the output matrix
* @param mat1 Pointer to first element of the 1st 2D matrix
* @param mat2 Pointer to first element of the 2nd 2D matrix
* @param size Size of the matrix
*/
double* mult_mat2D_brute(double* mat1, double* mat2, int row1, int col1, int row2, int col2);

double* tr033(double* mat, int rows, int cols);
const double* tr033(const double* mat, int rows, int cols);

/*
* Element-wise multiplies two 2D matricies with the same size
* @param mat_out Pointer to first element of the output matrix
* @param mat1 Pointer to first element of the 1st 2D matrix
* @param mat2 Pointer to first element of the 2nd 2D matrix
* @param size Size of the matrix
*/
void emult_mat2D(double* p_mat_out, double* p_mat1, double* p_mat2, int size);

/*
* Element-wise multiplies two 2D matricies with the same size. Returns a pointer to a 1D matrix
* @param mat1 Pointer to first element of the 1st 2D matrix
* @param mat2 Pointer to first element of the 2nd 2D matrix
* @param size Size of the matrix
* @param row Number of rows
*/
double** emult_mat2D(double* mat1, double* mat2, int size, int row);

/*
* Scales a 2D matrix by some value 
* @param mat Pointer to first element of the 2D matrix
* @param size Size of the matrix
* @param val Scalar
*/
void scale_mat2D(double* p_mat, int size, double val);

/*
* Outputs the data stored in a 2D array to the terminal
* @param mat Pointer to the first element of the 2D matrix
* @param col Number of columns
* @param row Number of rows
*/
void out_mat2D(double* p_mat, int row, int col);

/*
* Outputs the data stored in a 2D array to the terminal
* @param mat Pointer to the first element of the 2D matrix
* @param col Number of columns
* @param row Number of rows
*/
void out_mat2D(const double* p_mat, int row, int col);

void fill_mat2D(double* p_mat, int row, int col, double val);

double* p_mult_mat2D(double* p_mat1, double* p_mat2, int row1, int col1, int page, int row2, int col2);

double* p_emult_mat2D(double* p_mat1, double* p_mat2, int page, int row, int col);

double* p_cmult_mat2D(double* p_mat, int page, int row, int col, double val);

double* p_cadd_mat2D(double* p_mat, int page, int row, int col, double val);

double* p_sum_mat2D(double* p_mat, int page, int row, int col);

double* p_emult_mat2D_a(double* p_mat1, double* p_mat2, int page, int row, int col);

double* p_eadd_mat2D(double* p_mat1, double* p_mat2, int page, int row, int col);

double* p_eadd_mat2D_b(double* p_mat1, double* p_mat2, int page, int row, int col);

double* p_emult_mat2D_b(double* p_mat1, double* p_mat2, int page, int row, int col);

#endif