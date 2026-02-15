#include <iostream>
#include <vector>

#include <grid.cpp>
#include <field.cpp>

using namespace std;

void apply_stencil(const Grid &g,const Field &u, Field &Lu){
    int Nx=g.Nx,Ny=g.Ny;
    double hx2 = g.hx*g.hx;
    double hy2 = g.hy*g.hy;

    for(int j=1;j<=Ny;j++){
    for(int i=1;i<=Nx;i++){
      double uC = u.atc(i,j);
      double ux = (u.atc(i+1,j) - 2.0*uC + u.atc(i-1,j)) / hx2;
      double uy = (u.atc(i,j+1) - 2.0*uC + u.atc(i,j-1)) / hy2;
      Lu.at(i,j) = ux + uy;
    }
  }
}