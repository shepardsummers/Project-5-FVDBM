#include <iostream>
#include <iomanip>
#include <cmath>
#include "smm.h"
#include <fstream>

using namespace std;

void startup();
void write_mat(const string& filename, double* mat, int rows, int cols);
void write_mat(const string& filename, const double* mat, int rows, int cols);
void write_mat(const string& filename, double* mat, int pages, int rows, int cols);
double* eqm_d2q9(double* p_Rho, double* p_U, double* p_ksi, double* p_w, int row, int col);

int main() {
    startup();
    cout << "Initializing variables" << endl;
    // Node matrix size
    const int N_nd_x = 100;
    const int N_nd_y = N_nd_x;
    // Cell matrix size;
    const int N_cv_x = N_nd_x - 1;
    const int N_cv_y = N_nd_y - 1;
    // Total size
    const int L = N_nd_x;
    const int H = N_nd_y;
    // Dist between nodes
    const int dx = L/N_nd_x;
    const int dy = H/N_nd_y;
    const double A = dx*dy; // Cell area
    const int cv_length[4] = {dy, dx, dy, dx}; // Cell dimensions
    // Timestep
    const double dt = 0.1;
    // Lattice velocity
    double ksi[18] = {0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,-1,-1,1,-1};
    double* p_ksi = &ksi[0];
    int ksi_row = 9;
    int ksi_col = 2;
    write_mat("output/ksi.dat", &ksi[0], 9, 2);
    // Weights
    double w[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    double* p_w = &w[0];
    // Fluid info
    const int Re = 100; // Reanolds number
    const double tau = 0.4; // Relaxation time
    const double c_s = 1/sqrt(3); // Speed of sound for D2Q9
    const double vis = pow(tau*c_s, 2); // Kinematic viscosity
    const double U_top[2] = {{Re*vis/L},{0}}; // Top velocity
    cout << "Top velocity: " << setprecision(15) << U_top[0] << endl;
    const double rho_init = 2; // Initial density throughout volume
    // Initializing arrays
    const int f_emp[9] = {0}; // Empty pdf
    double Rho_nd[N_nd_y][N_nd_x]; // Density at nodes
    fill_mat2D(&Rho_nd[0][0], N_nd_y, N_nd_x, rho_init);
    double Rho_cv[N_cv_y][N_cv_x]; // Density in cells
    fill_mat2D(&Rho_cv[0][0], N_cv_y, N_cv_x, rho_init);
    double U_nd[2][N_nd_y][N_nd_x] = {0}; // Velocity at nodes
    double U_cv[2][N_cv_y][N_cv_x] = {0}; // Velocity in cells

    double U_test[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
    double* p_U_test = &U_test[0];
    int U_test_page = 2;
    int U_test_row = 2;
    int U_test_col = 3;
    write_mat("output/U_test.dat", p_U_test, U_test_page, U_test_row, U_test_col);

    double* p_ksiU = p_mult_mat2D(p_ksi, p_U_test, ksi_row, ksi_col, U_test_page, U_test_row, U_test_col);
    int ksiU_page = ksi_row;
    int ksiU_row = U_test_row;
    int ksiU_col = U_test_col;
    write_mat("output/ksiU_test.dat", p_ksiU, ksiU_page, ksiU_row, ksiU_col);

    double* ksiU_2 = p_emult_mat2D(p_ksiU, p_ksiU, 9, 2, 3);
    write_mat("output/ksiU_2_test.dat", ksiU_2, 9, 2, 3);

    double* erm = p_sum_mat2D(p_U_test, U_test_page, U_test_row, U_test_col);
    write_mat("output/U_sum.dat", erm, 2, 3);

    double* test = p_emult_mat2D_a(p_ksiU, p_w, 9, 2,3);
    write_mat("output/test.dat", test, 9, 2, 3);

    double rho_test[6] = {2, 2, 2, 2, 2, 2};

    double* out = eqm_d2q9(rho_test, p_U_test, p_ksi, p_w, 2, 3);

    return 0;
}

void startup() {
    cout << "+----------------------------------+" << endl
         << "FVDBM LDC SOLVER V0.1" << endl
         << "Data output in .dat file" << endl
         << "+----------------------------------+" << endl;
}

void write_mat(const string& filename, double* mat, int rows, int cols) {
    ofstream out(filename);
    if (!out) {
        cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            // << fixed <-- will fix it to do decimals and not scientific notation
            out  << setprecision(15) << mat[j*cols + i];
            if (i < cols - 1)
                out << " ";  // Space-separated values
        }
        out << "\n";  // Newline at the end of each row
    }

    out.close();
}

void write_mat(const string& filename, const double* mat, int rows, int cols) {
    ofstream out(filename);
    if (!out) {
        cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            // << fixed <-- will fix it to do decimals and not scientific notation
            out  << setprecision(15) << mat[j*cols + i];
            if (i < cols - 1)
                out << " ";  // Space-separated values
        }
        out << "\n";  // Newline at the end of each row
    }

    out.close();
}

void write_mat(const string& filename, double* mat, int pages, int rows, int cols) {
    ofstream out(filename);
    if (!out) {
        cerr << "Failed to open file: " << filename << "\n";
        return;
    }
    for (int p = 0; p < pages; ++p) {
        out << "page " << p+1 << "\n";
        for (int j = 0; j < rows; ++j) {
            for (int i = 0; i < cols; ++i) {
                // << fixed <-- will fix it to do decimals and not scientific notation
                out  << setprecision(15) << mat[i + (p*rows + j)*cols];
                if (i < cols - 1)
                    out << " ";  // Space-separated values
            }
            out << "\n";  // Newline at the end of each row
        }
    }

    out.close();
}

double* eqm_d2q9(double* p_Rho, double* p_U, double* p_ksi, double* p_w, int row, int col) {
    
    double* p_ksiU = p_mult_mat2D(p_ksi, p_U, 9, 2, 2, row, col);
    double* p_out = p_cmult_mat2D(p_ksiU, 9, row, col, 3.0);
    p_out = p_cadd_mat2D(p_out, 9, row, col, 1.0);
    double* p_out2 = p_emult_mat2D(p_ksiU, p_ksiU, 9, row, col);
    p_out2 = p_cmult_mat2D(p_out2, 9, row, col, (9.0/2.0));
    p_out = p_eadd_mat2D(p_out, p_out2, 9, row, col);
    double* p_out3 = p_sum_mat2D(p_emult_mat2D(p_U, p_U, 2, row, col), 2, row, col);
    p_out3 = p_cmult_mat2D(p_out3, 1, row, col, (-3.0/2.0));
    
    p_out = p_eadd_mat2D_b(p_out, p_out3, 9, row, col);
    p_out = p_emult_mat2D_a(p_emult_mat2D_b(p_out, p_Rho, 9, row, col), p_w, 9, row, col);

    return p_out;
};