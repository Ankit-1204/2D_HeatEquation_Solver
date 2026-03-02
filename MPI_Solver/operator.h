#pragma once
#include <cassert>
#include "grid.h"
#include "field.h"

inline void apply_A_local(const LocalGrid &lg, const Field &x, Field &Ax) {
    assert(x.g == Ax.g && x.g == &lg);
    for (int j = 1; j <= lg.Ny_loc; ++j) {
        for (int i = 1; i <= lg.Nx_loc; ++i) {
            double c = x.atc(i,j);
            Ax.at(i,j) = 4.0*c - x.atc(i+1,j) - x.atc(i-1,j) - x.atc(i,j+1) - x.atc(i,j-1);
        }
    }
}

inline void apply_M_local(const LocalGrid &lg, const Field &x, Field &y, double theta, double gamma) {
    Field Ax(&lg); 
    Ax.reset(&lg); 
    apply_A_local(lg, x, Ax);
    for (int j = 1; j <= lg.Ny_loc; ++j)
        for (int i = 1; i <= lg.Nx_loc; ++i)
            y.at(i,j) = x.atc(i,j) + theta * gamma * Ax.atc(i,j);
}