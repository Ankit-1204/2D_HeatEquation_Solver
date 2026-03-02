#pragma once
#include <cassert>


struct LocalGrid {
    int Nx_loc = 0, Ny_loc = 0; 
    int gx = 0, gy = 0;         
    double x0 = 0.0, x1 = 1.0, y0 = 0.0, y1 = 1.0;
    double hx = 0.0, hy = 0.0; 

    LocalGrid() = default;
    LocalGrid(double x0_, double x1_, double y0_, double y1_,
              int Nx_loc_, int Ny_loc_, int gx_, int gy_) {
        init(x0_, x1_, y0_, y1_, Nx_loc_, Ny_loc_, gx_, gy_);
    }
    void init(double x0_, double x1_, double y0_, double y1_,
              int Nx_loc_, int Ny_loc_, int gx_, int gy_) {
        x0 = x0_; x1 = x1_; y0 = y0_; y1 = y1_;
        Nx_loc = Nx_loc_; Ny_loc = Ny_loc_;
        gx = gx_; gy = gy_;
    }

    inline int stride() const { return Nx_loc + 2; } 
    inline int idx(int i,int j) const { return i + j * stride(); } 

    
    inline double global_x(int i_local) const {
        assert(i_local >= 1 && i_local <= Nx_loc);
        int global_interior_index = gx + i_local; 
        return x0 + global_interior_index * hx;
    }
    inline double global_y(int j_local) const {
        assert(j_local >= 1 && j_local <= Ny_loc);
        int global_interior_index = gy + j_local;
        return y0 + global_interior_index * hy;
    }
};