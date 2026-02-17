#include <iostream>
#include <vector>

#include <grid.cpp>
#include <field.cpp>
#include <grid.h>
#include <field.h>

using namespace std;

inline void for_interior(const Grid &g, function<void(int,int)> fn) {
    for (int j = 1; j <= g.Ny; ++j)
        for (int i = 1; i <= g.Nx; ++i)
            fn(i,j);
}

inline double dot_interior(const Grid &g, const Field &a, const Field &b) {
    double s = 0.0;
    for_interior(g, [&](int i,int j){ s += a.atc(i,j) * b.atc(i,j); });
    return s;
}

inline void copy_interior(const Grid &g, const Field &src, Field &dst) {
    for_interior(g, [&](int i,int j){ dst.at(i,j) = src.atc(i,j); });
}
inline void axpy_interior(const Grid &g, double alpha, const Field &x, Field &y) {
    for_interior(g, [&](int i,int j){ y.at(i,j) += alpha * x.atc(i,j); });
}
void apply_stencilA(const Grid &g,const Field &u, Field &Au){
    int Nx=g.Nx,Ny=g.Ny;
    double hx2 = g.hx*g.hx;
    double hy2 = g.hy*g.hy;

    for(int j=1;j<=Ny;j++){
    for(int i=1;i<=Nx;i++){
      double uC = u.atc(i,j);
      double ux = (u.atc(i+1,j) - 2.0*uC + u.atc(i-1,j)) / hx2;
      double uy = (u.atc(i,j+1) - 2.0*uC + u.atc(i,j-1)) / hy2;
      Au.at(i,j) = ux + uy;
    }
  }
}

void apply_M(const Grid &g, const Field &u,Field &y,double theta, double gamma){
  Field Au(g);
  apply_stencilA(g,u,Au);
  int N=g.Nx*g.Ny;
  for(int i=1;i<=g.Nx;i++){
    for(int j=1;j<=g.Ny;j++){
      y.at(i,j)=u.atc(i,j) + theta*gamma*Au.atc(i,j);
    }
  }
}

inline void fill_ghosts_dirichlet(Field &u,
                                  function<double(double)> left, 
                                  function<double(double)> right, 
                                  function<double(double)> bottom,
                                  function<double(double)> top)   
{
    const Grid &g = u._grid;
    
    for (int j=1;j<=g.Ny;++j) {
        double y = g.y(j);
        u.at(0,j) = left(y);
        u.at(g.Nx+1,j) = right(y);
    }

    for (int i=0;i<=g.Nx+1;++i) {
        double x = g.x(i);
        u.at(i,0) = bottom(x);
        u.at(i,g.Ny+1) = top(x);
    }
}
inline double default_source(double x, double y, double t) { return 0.0; }

inline void build_rhs_BE(const Grid &g, const Field &u_old, Field &b, double dt,
                         function<double(double,double,double)> source, double t_next)
{
    for_interior(g, [&](int i,int j){
        double x = g.x(i); double y = g.y(j);
        b.at(i,j) = u_old.atc(i,j) + dt * source(x,y,t_next);
    });
}