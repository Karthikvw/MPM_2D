#pragma once 

#include<iostream>
#include<array>
#include<vector>
#include<cmath>
#include<functional>
#include<cfloat>
#include<utility>

#include<output_util.hpp>
#include<MPM_Grid.hpp>
#include<MPM_Body.hpp>

using namespace std;

class MPM_Solver
{
    public:
    //Constructor
    MPM_Solver(){}

    //Solution function
    void Solve(MPM_Grid* Grid, vector<MPM_Body*> Bodies, size_t NoS, double dt);

    //Boundary conditions
    public:
    class BC
    {
        public:
        BC(const function<bool(array<double,2>)> &Criterion_) : Criterion(Criterion_), active_vx(false), active_vy(false) {}
        
        const function<bool(array<double,2>)> Criterion;
        bool   active_vx;       double value_vx;
        bool   active_vy;       double value_vy;
    };

    //Container of the solver class to carry boundary conditions
    vector<BC> BoundaryConditions;
    
    public:
    //Destructor
    ~MPM_Solver(){}
};

