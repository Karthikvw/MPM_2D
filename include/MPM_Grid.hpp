#pragma once 

#include<iostream>
#include<array>
#include<vector>
#include<string>
#include<cmath>
#include<output_util.hpp>
#include<para_export.hpp>

using namespace std;

class MPM_Grid
{
    public:
    vector<array<double, 11>> NC;   //Nodal container  
    double lx,ly,dx,dy,x_0,y_0;     //length of gird, size of cells and origin of the grid
    size_t nx,ny,NoNodes;           //number of cells

    public:
    //constructor
    MPM_Grid(array<double,2> L, array<size_t,2> N,array<double,2> O={0.0,0.0})
    /*    
        Initialises the 2D grid

        Example: MPMGrid_2D([4,4], [5,5], [1,1])
    
        The first argument is length of grid in X and Y direction
    
        The second argument is number of cells in X and Y direction
    
        The third argument is Origin and is optional

        *** Default Origin - Spatial Origin [0,0]***
    */

    {
        lx = L[0]; ly = L[1];
        nx = N[0]; ny = N[1];
        dx = lx/nx; dy = ly/ny;
        x_0 = O[0]; y_0 = O[1];
        NoNodes = (nx+1)*(ny+1);
        NC.resize(NoNodes);
    }
    
    public:
    array<array<double, 2>, 4> Cell_X(array<double,2> X);

    public:
    array<size_t,4> Cell_I(array<double,2> X);
     
    public:
    array<double,12> SHP(array<double,2> X);

    public:
    void ResetNC();

    public:
    void SolveGrid(double dt, double alpha, double mass_cutoff = 1e-10 );

    public:
    vector<array<double, 2>> Grid_X();

    public:
    void VTKGrid(string output_file_name, size_t step);

    public:
    //destructor
    ~MPM_Grid(){}

};