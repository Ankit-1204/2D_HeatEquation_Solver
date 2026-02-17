#include <iostream>
#include <vector>

#include <grid.h>
using namespace std;

class Field{
    public:
        Grid _grid;
        vector<double> fields;
    Field(const Grid& grid) :_grid(grid) {
        fields.resize((grid.Nx+2)*(grid.Ny+2),0);
    }
    inline double& at(int i,int j)            { return fields[_grid.idx(i,j)]; } // (A)
    inline double  atc(int i,int j) const    { return fields[_grid.idx(i,j)]; } // (B)

};

