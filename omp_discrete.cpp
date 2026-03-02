#include <iostream>
#include "field.h"
#include "grid.h"

void apply_A(const Grid &g,const Field &x, Field &Ax){
        #pragma omp parallel for schedule(static)
         for (int j=1;j<=g.Ny;++j) {
        int stride = g.stride;
        int base = j*stride;
        for (int i=1;i<=g.Nx;++i) {
            double c = x.fields[base + i];
            double sum = 4.0*c - x.fields[base + (i+1)] - x.fields[base + (i-1)]
                                - x.fields[(j+1)*stride + i] - x.fields[(j-1)*stride + i];
            Ax.fields [base + i] = sum;
        }
    }
}
inline void apply_M(const Grid &g, const Field &x, Field &y, double theta, double gamma) {
    Field Ax(g); 
    apply_A(g, x, Ax);
    #pragma omp parallel for schedule(static)
    for (int j=1;j<=g.Ny;++j) {
        int stride = g.stride;
        int base = j*stride;
        for (int i=1;i<=g.Nx;++i)
            y.fields[base + i] = x.fields[base + i] + theta * gamma * Ax.fields[base + i];
    }
}

inline double dot(const Grid &g, const Field &a, const Field &b) {
    double s = 0.0;
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (int j=1;j<=g.Ny;++j) {
        int stride = g.stride;
        int base = j*stride;
        for (int i=1;i<=g.Nx;++i) s += a.fields[base + i] * b.fields[base + i];
    }
    return s;
}

inline void axpy(const Grid &g, double alpha, const Field &x, Field &y) {
    #pragma omp parallel for schedule(static)
    for (int j=1;j<=g.Ny;++j) {
        int stride = g.stride;
        int base = j*stride;
        for (int i=1;i<=g.Nx;++i) y.fields[base + i] += alpha * x.fields[base + i];
    }
}