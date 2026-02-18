#pragma once
#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include <functional>
#include "grid.h"
#include "field.h"

struct CGResult {
    int iters;
    double resnorm;
};

CGResult cg_solver(
    const Grid &g,
    std::function<void(const Field&, Field&)> applyM,
    const Field &b,
    Field &x,
    std::function<void(Field&)> fill_ghosts_for_x,
    int maxit = 1000,
    double tol = 1e-8
);


#endif 
