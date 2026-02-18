#include "field.h"

Field::Field(const Grid& grid) : _grid(grid) {
    fields.resize((grid.Nx + 2) * (grid.Ny + 2), 0.0);
}

double& Field::at(int i, int j) {
    return fields[_grid.idx(i, j)];
}

double Field::atc(int i, int j) const {
    return fields[_grid.idx(i, j)];
}
