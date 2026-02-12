#include <iostream>

class Grid{
    public :
      int Nx,Ny;
      double xmin,xmax,ymin,ymax,hx,hy;

    Grid(double xmin_, double xmax_, double ymin_, double ymax_, int Nx_, int Ny_)
   : Nx(Nx_), Ny(Ny_), xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_) {
     hx = (xmax-xmin)/(Nx+1);
     hy = (ymax-ymin)/(Ny+1);
  }
};