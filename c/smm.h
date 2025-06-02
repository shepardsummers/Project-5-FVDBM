#ifndef SMM_H
#define SMM_H

/*
* Scales a 2D matrix by some value 
* @param mat Pointer to first element of the 2D matrix
* @param size Size of the matrix
* @param val Scalar
*/
void scale_mat2D(double* mat, int size, double val);

/*
* Element-wise multiplies two 2D matricies with the same size
* @param mat_out Pointer to first element of the output matrix
* @param mat1 Pointer to first element of the 1st 2D matrix
* @param mat2 Pointer to first element of the 2nd 2D matrix
* @param size Size of the matrix
*/
void emult_mat2D(double* mat_out, double* mat1, double* mat2, int size);

/*
* Outputs the data stored in a 2D array to the terminal
* @param mat Pointer to the first element of the 2D matrix
* @param col Number of columns
* @param row Number of rows
*/
void out_mat2D(double* mat, int col, int row);

/*
* Outputs the data stored in a 2D array to the terminal
* @param mat Pointer to the first element of the 2D matrix
* @param col Number of columns
* @param row Number of rows
*/
void out_mat2D(const double* mat, int col, int row);

#endif