#pragma once
#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

#include <functional>  
#include "grid.h"
#include "field.h"

void for_interior(const Grid &g, std::function<void(int,int)> fn);

double dot_interior(const Grid &g, const Field &a, const Field &b);

void copy_interior(const Grid &g, const Field &src, Field &dst);

void axpy_interior(const Grid &g, double alpha, const Field &x, Field &y);

void apply_stencilA(const Grid &g, const Field &u, Field &Au);

void apply_M(const Grid &g, const Field &u, Field &y, double theta, double gamma);

void fill_ghosts_dirichlet(Field &u,
                           std::function<double(double)> left,
                           std::function<double(double)> right,
                           std::function<double(double)> bottom,
                           std::function<double(double)> top);

double default_source(double x, double y, double t);

void build_rhs_BE(const Grid &g, const Field &u_old, Field &b, double dt,
                  std::function<double(double,double,double)> source, double t_next);

#endif
