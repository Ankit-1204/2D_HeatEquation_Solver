#include <iostream>
#include "field.h"
#include "grid.h"
#include "discretization.h"
#include "cg_solver.h"
#include <cmath>

using namespace std;

void main(int argc, char** argv){
    int Nx = 50, Ny = 50;      
    double x0=0.0, x1=1.0, y0=0.0, y1=1.0;
    double alpha = 1.0;       
    double dt = 0.0005;            
    double t_final = 0.01;
    double theta = 0.5;            
    int output_every = 10;
    double tol = 1e-8;
    int maxit = 2000;

    Grid g(x0,x1,y0,y1,Nx,Ny);
    Field u_old(g), u_new(g), b(g), temp(g);

    // for easy dirchlet Bc
    auto g_left = [&](double y){ 
        return 0.0;};
    auto g_right = [&](double y){
        return 0.0;};
    auto g_bottom= [&](double x){
        return 0.0;};
    auto g_top = [&](double x){
        return 0.0;};

    // ezy initial condition
    double cx = 0.5*(x0+x1), cy = 0.5*(y0+y1);
    double sigma = 0.08;
    for_interior(g, [&](int i,int j){
        double x = g.x(i), y = g.y(j);
        u_old.at(i,j) = exp(-((x-cx)*(x-cx)+(y-cy)*(y-cy))/(2.0*sigma*sigma));
    });
    fill_ghosts_dirichlet(u_old, g_left, g_right, g_bottom, g_top);
    
    double hx = g.hx, hy = g.hy;
    double h = hx; 
    double gamma = alpha * dt / (h*h);
     auto applyM_lambda = [&](const Field &x, Field &y) {
        apply_M(g, x, y, theta, gamma);
    };
    auto fill_ghosts_for_x = [&](Field &xfield){
        fill_ghosts_dirichlet(xfield, g_left, g_right, g_bottom, g_top);
    };
    int step = 0; double t = 0.0;
    int nsteps = (int)ceil(t_final / dt);

    for (step = 0; step < nsteps; ++step) {
        double t_next = t + dt;

        fill_ghosts_dirichlet(u_old, g_left, g_right, g_bottom, g_top);
        build_rhs_BE(g, u_old, b, dt, default_source, t_next);

        copy_interior(g, u_old, u_new);
        fill_ghosts_for_x(u_new);
        auto applyM_fn = [&](const Field &xin, Field &yout){
            applyM_lambda(xin, yout);
        };

        CGResult res = cg_solver(g, applyM_fn, b, u_new,
                                 [&](Field &xf){ fill_ghosts_dirichlet(xf, g_left, g_right, g_bottom, g_top); },
                                 maxit, tol);
        
        if (step % output_every == 0) {
            cout << "step="<<step<<" t="<<t_next<<" pcg_iters="<<res.iters<<" res="<<res.resnorm<<"\n";       
        }
        u_old = u_new;
        t=t_next;
    }
    cout << "Done. final time t="<<t<<"\n";
    return ;
}
