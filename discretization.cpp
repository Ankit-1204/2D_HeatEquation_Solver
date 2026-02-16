#include <iostream>
#include <vector>

#include <grid.cpp>
#include <field.cpp>

using namespace std;

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

void apply_M(const Grid &g, const Field &u,Field &Au,Field &y,double theta, double gamma){
  apply_stencilA(g,u,Au);
  int N=g.Nx*g.Ny;
  for(int i=1;i<=g.Nx;i++){
    for(int j=1;j<=g.Ny;j++){
      y.at(i,j)=u.atc(i,j) + theta*gamma*Au.atc(i,j);
    }
  }
}