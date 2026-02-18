#include <iostream>
#include "field.h"
#include "grid.h"

void apply_A(const Grid &g,const Field &x, Field &Ax){
        #pragma omp parallel for schedule(static)
         for (int j=1;j<=g.Ny;++j) {
        int stride = g.stride();
        int base = j*stride;
        for (int i=1;i<=g.Nx;++i) {
            double c = x.field[base + i];
            double sum = 4.0*c - x.field[base + (i+1)] - x.field[base + (i-1)]
                                - x.field[(j+1)*stride + i] - x.field[(j-1)*stride + i];
            Ax.field [base + i] = sum;
        }
    }
}
inline void apply_M(const Grid &g, const Field &x, Field &y, double theta, double gamma) {
    Field Ax(&g); 
    apply_A(g, x, Ax);
    #pragma omp parallel for schedule(static)
    for (int j=1;j<=g.Ny;++j) {
        int stride = g.stride();
        int base = j*stride;
        for (int i=1;i<=g.Nx;++i)
            y.d[base + i] = x.d[base + i] + theta * gamma * Ax.d[base + i];
    }
}