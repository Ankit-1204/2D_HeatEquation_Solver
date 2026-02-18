#include "grid.h"

Grid::Grid(double xmin_, double xmax_,
           double ymin_, double ymax_,
           int Nx_, int Ny_)
    : Nx(Nx_), Ny(Ny_),
      xmin(xmin_), xmax(xmax_),
      ymin(ymin_), ymax(ymax_)
{
    hx = (xmax - xmin) / (Nx + 1);
    hy = (ymax - ymin) / (Ny + 1);
    stride = Nx + 2;
}

int Grid::idx(int i, int j) const {
    return i + stride * j;
}

double Grid::x(int i) const {
    return xmin + i * hx;
}

double Grid::y(int j) const {
    return ymin + j * hy;
}
