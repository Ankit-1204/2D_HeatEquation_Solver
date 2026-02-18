#pragma once
#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include "grid.h"

class Field {
public:
    Grid _grid;
    std::vector<double> fields;

    explicit Field(const Grid& grid);

    double& at(int i, int j);     
    double  atc(int i, int j) const; 
};

#endif
