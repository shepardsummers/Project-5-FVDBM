#include <iostream>
#include <iomanip>
#include <cmath>
#include "smm.h"
#include <fstream>

using namespace std;

void startup();
void write_mat(const string& filename, double* mat, int rows, int cols);

int main() {
    startup();
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
    const double A = dx*dy;
    // Timestep
    const double dt = 0.1;
    // Lattice velocity
    const double ksi[2][9] = {{0, 1, 0,-1, 0, 1,-1,-1, 1},
                              {0, 0, 1, 0,-1, 1, 1,-1,-1}};
    // Weights
    const double w[9] = {4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36};
    // Fluid info
    const int Re = 100; // Reanolds number
    const double tau = 0.4; // Relaxation time
    const double c_s = 1/sqrt(3); // Speed of sound for D2Q9
    const double vis = pow(tau*c_s, 2); // Kinematic viscosity
    const double U_top[2] = {{Re*vis/L},{0}}; // Top velocity
    const double rho_init = 2; // Initial density throughout volume
    // Initializing arrays
    const int f_emp[9] = {0}; // Empty pdf
    double Rho_nd[N_nd_y][N_nd_x]; // Density at nodes
    //fill_n(Rho_nd, N_nd_y*N_nd_x, rho_init);
    double Rho_cv[N_cv_y][N_cv_x]; // Density in cells
    //fill_n(Rho_cv, N_cv_y*N_cv_x, rho_init);
    double U_nd[2][N_nd_y][N_nd_x] = {0}; // Velocity at nodes
    double U_cv[2][N_cv_y][N_cv_x] = {0}; // Velocity in cells


    cout << setprecision(15) << U_top[0] << endl;

    double matc[2][3];
    double mata[2][3] = {{1, 2, 3}, {4, 5, 6}};
    double matb[2][3] = {{1, 2, 3}, {4, 5, 6}};
    double matd[2][5] = {{1, 2, 3, 4, 5},{1, 2, 3, 4, 5}};
    double mate[5][2] = {{1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}};

    //emult_mat2D(&matc[0][0], &mata[0][0], &matb[0][0], 2*3);
    //double** out = emult_mat2D(&mata[0][0], &matb[0][0], 2*3, 3);

    double* guh = mult_mat2D_brute(&matd[0][0], &mate[0][0], 2, 5, 5, 2);

    //out_mat2D(guh, 2, 2);
    //out_mat2D(tr033(&mate[0][0], 5, 2), 2, 5);
    //out_mat2D(&ksi[0][0], 2, 9);
    write_mat("res.dat", guh, 2, 2);
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

    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < cols; i++) {
            out << mat[j * cols + i];
            if (i < cols - 1)
                out << " ";  // Space-separated values
        }
        out << "\n";  // Newline at the end of each row
    }

    out.close();
}