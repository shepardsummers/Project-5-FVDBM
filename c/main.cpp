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
    const double ksi[2][9] = {{0, 1, 0,-1, 0, 1,-1,-1, 1},
                              {0, 0, 1, 0,-1, 1, 1,-1,-1}};
    const double* ksi_T = tr033(&ksi[0][0], 2, 9);
    // Weights
    const double w[9] = {4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36};
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

    double guh[2][2][3] = {{{1, 2, 1}, {2, 1, 2}},
                           {{3, 4, 3}, {4, 0, 0}}};
    write_mat("guh.dat", &guh[0][0][0], 2, 2, 3);

    double erm[5][2] = {{1,6},{2,7},{3,8},{4,9},{5,10}};
    write_mat("erm.dat", &erm[0][0], 5, 2);

    double* super_guh = p_mult_mat2D(&erm[0][0], &guh[0][0][0], 5, 2, 2, 2, 3);
    write_mat("super_guh.dat", super_guh, 5, 2, 3);

    //write_mat("ksi_T.dat", ksi_T, 9, 2);
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

    for (int j = 0; j < rows; ++j) {
        out << "row " << j+1 << "\n";
        for (int i = 0; i < cols; ++i) {
            for (int p = 0; p < pages; ++p) {
                // << fixed <-- will fix it to do decimals and not scientific notation
                out  << setprecision(15) << mat[(j*cols + i)*pages + p];
                if (p < pages - 1)
                    out << " ";  // Space-separated values
            }
            out << "\n";  // Newline at the end of each col
        }
        out << "\n";  // Newline at the end of each row
    }

    out.close();
}