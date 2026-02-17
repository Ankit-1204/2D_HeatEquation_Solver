#include <iostream>
#include <grid.cpp>
#include <discretization.cpp>
#include <functional>
#include <cmath>

using namespace std;
struct CGResult { int iters; double resnorm; };
CGResult cg_solver(const Grid &g,function<void(const Field&, Field&)> applyM,const Field &b, Field &x,function<void(Field&)> fill_ghosts_for_x,
                   int maxit=1000, double tol=1e-8){
        Field r(g), z(g), p(g), Ap(g), temp(g);
        fill_ghosts_for_x(x);
        applyM(x, temp); 
        for(int i=1;i<=g.Nx;i++){
            for(int j=1;j<=g.Ny;j++){
            r.at(i,j) = b.atc(i,j) - temp.atc(i,j);
            }
        }
        double normb = sqrt(dot_interior(g,b,b));
        if (normb == 0.0) normb = 1.0;
        copy_interior(g, r, z);
        copy_interior(g, z, p);

        double rz_old = dot_interior(g, r, z);
        double res0 = sqrt(rz_old);
        double resnorm = res0;
        if (res0 <= tol * normb) return {0, res0};
        int iter=0;
        for(iter =1;iter<=maxit;iter++){
            fill_ghosts_for_x(p);
            applyM(p, Ap);
            double pAp = dot_interior(g, p, Ap);
            if (fabs(pAp) < 1e-18) break;
            double alpha = rz_old / pAp;

            axpy_interior(g, alpha, p, x);
            fill_ghosts_for_x(x);
            axpy_interior(g, -alpha, Ap, r);

            double res = sqrt(dot_interior(g,r,r));;
            if (res <= tol * normb) {
                resnorm = res; break;
            }
            copy_interior(g, r, z);
            double rz_new = dot_interior(g, r, z);
            double beta = rz_new / rz_old;
            for_interior(g, [&](int i,int j){ p.at(i,j) = z.atc(i,j) + beta * p.atc(i,j); });
            rz_old = rz_new;
        }
        return {iter, resnorm};
}

inline void build_rhs_BE(const Grid &g, const Field &u_old, Field &b, double dt,
                         function<double(double,double,double)> source, double t_next)
{
    for_interior(g, [&](int i,int j){
        double x = g.x(i); double y = g.y(j);
        b.at(i,j) = u_old.atc(i,j) + dt * source(x,y,t_next);
    });
}