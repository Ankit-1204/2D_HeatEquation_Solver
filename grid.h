#pragma once
#ifndef GRID_H
#define GRID_H

class Grid {
public:
    int Nx, Ny;
    double xmin, xmax, ymin, ymax;
    double hx, hy;
    int stride;

    Grid(double xmin_, double xmax_,
         double ymin_, double ymax_,
         int Nx_, int Ny_);

    int idx(int i, int j) const;

    double x(int i) const;  
    double y(int j) const;
};

#endif
