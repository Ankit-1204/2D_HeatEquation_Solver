#include <iostream>
#include <vector>
#include <mpi.h>
#include <functional>
using namespace std;

inline double bc_left(double y)   { return 0.0; }
inline double bc_right(double y)  { return 0.0; }
inline double bc_bottom(double x) { return 0.0; }
inline double bc_top(double x)    { return 0.0; }

inline double global_x(int global_i, double x0, double hx) { return x0 + (global_i + 1) * hx; }
inline double global_y(int global_j, double y0, double hy) { return y0 + (global_j + 1) * hy; }

inline void apply_A_local(const LocalGrid &lg, const Field &x, Field &Ax) {
    for_interior_loc(lg, [&](int i,int j){
        double c = x.atc(i,j);
        double sum = 4.0*c - x.atc(i+1,j) - x.atc(i-1,j) - x.atc(i,j+1) - x.atc(i,j-1);
        Ax.at(i,j) = sum;
    });
}

inline void for_interior_loc(const LocalGrid &g, function<void(int,int)> fn) {
    for (int j = 1; j <= g.Ny_loc; ++j)
        for (int i = 1; i <= g.Nx_loc; ++i)
            fn(i,j);
}