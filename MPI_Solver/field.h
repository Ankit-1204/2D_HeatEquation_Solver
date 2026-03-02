#pragma once
#include <vector>
#include "grid.h"

using namespace std;
struct Field {
    const LocalGrid *g = nullptr;
    vector<double> d; 
    Field() = default;
    Field(const LocalGrid *grid) { reset(grid); }
    void reset(const LocalGrid *grid) {
        g = grid;
        d.assign((g->Nx_loc + 2) * (g->Ny_loc + 2), 0.0);
    }
    inline double& at(int i,int j) { return d[g->idx(i,j)]; }
    inline double  atc(int i,int j) const { return d[g->idx(i,j)]; }
};