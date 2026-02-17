#include <iostream>

class Grid{
    public :
      int Nx,Ny;
      double xmin,xmax,ymin,ymax,hx,hy;
      int stride;

    Grid(double xmin_, double xmax_, double ymin_, double ymax_, int Nx_, int Ny_)
   : Nx(Nx_), Ny(Ny_), xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_) {
     hx = (xmax-xmin)/(Nx+1);
     hy = (ymax-ymin)/(Ny+1);
     stride=Nx+2;
  }

  inline int idx(int i,int j)const {
    return i+ (stride*j);
  }
  inline double x(int i) const { return xmin + i * hx; } // note: i=0..Nx+1 where i=0 is ghost at x0
  inline double y(int j) const { return ymin + j * hy; }
};

