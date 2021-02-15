#include <MPM_Grid.hpp>

array<array<double, 2>, 4> MPM_Grid::Cell_X(array<double,2> X)    
/*    
    Returns local nodal coordinates of the given material point's coordinates

    Syntax: Cell_X(X)

    X - Coordinates of material point
        
    Example: Cell_X([1,3])
*/
{
    array<array<double, 2>, 4> X_I = {{{0,0}, {0,0}, {0,0}, {0,0}}};

    if (X[0] < x_0 || X[1] < y_0)
        throw runtime_error("Material point outside the grid");

    else if(X[0] > x_0+lx || X[1] > y_0+ly)
        throw runtime_error("Material point outside the grid");

    else if(X[0] == x_0+lx || X[1] == y_0+ly)
    {
        X_I[0][0] = (floor((X[0]-x_0)/dx) * dx) + x_0;
        X_I[0][1] = (floor((X[1]-y_0)/dy) * dy) + y_0;

        X_I[1][0] = X_I[0][0] - dx;
        X_I[1][1] = X_I[0][1];
            
        X_I[2][0] = X_I[0][0] - dx;
        X_I[2][1] = X_I[0][1] - dy;

        X_I[3][0] = X_I[0][0];
        X_I[3][1] = X_I[0][1] - dy;

        return X_I;
    }

    else
    {
        X_I[0][0] = (floor((X[0]-x_0)/dx) * dx) + x_0;
        X_I[0][1] = (floor((X[1]-y_0)/dy) * dy) + y_0;

        X_I[1][0] = X_I[0][0] + dx;
        X_I[1][1] = X_I[0][1];
            
        X_I[2][0] = X_I[0][0] + dx;
        X_I[2][1] = X_I[0][1] + dy;

        X_I[3][0] = X_I[0][0];
        X_I[3][1] = X_I[0][1] + dy;
        
        return X_I;
    }
}

array<size_t,4> MPM_Grid::Cell_I(array<double,2> X)
/*
    Returns local nodal IDs of the given material point's coordinates
        
    Syntax: Cell_X(X)

    X - Coordinates of material point
        
    Example: Cell_I([1,3])
*/
{
    array<size_t,4> N_I={0,0,0,0};

    if (X[0] < x_0 || X[1] < y_0)  
        throw runtime_error("Material point outside the grid");      

    else if(X[0] > x_0+lx || X[1] > y_0+ly)
        throw runtime_error("Material point outside the grid");

    else
    {
        size_t nxi,nyi;

        nxi = floor((X[0]-x_0)/dx);
        nyi = floor((X[1]-y_0)/dy);

        N_I[0] = nyi*(nx+1) + nxi;
        N_I[1] = nyi*(nx+1) + (nxi+1);
        N_I[2] = (nyi+1)*(nx+1) + (nxi+1);
        N_I[3] = (nyi+1)*(nx+1) + nxi;

        return N_I;
    }
}

array<double,12> MPM_Grid::SHP(array<double,2> X)
/*
    Returns a vector of shape functions and derived shape functions of the given material point (SHPC_mp)
        
    Syntax: SHP(X)

    X - Coordinates of material point
                
    Example: SHP([1,3])

    SHPC -> Stores Shape functions and derived shape functions of all material points in the given order

    Node 1 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |,

    Node 2 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |,

    Node 3 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |,

    Node 4 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |.
*/
{
    array<double,12> SHPC_mp;
    for (size_t i = 0; i < 12; i++)
    {
        SHPC_mp[i] = 0.0;
    }

    if (X[0] < x_0 || X[1] < y_0)
        throw runtime_error("Material point outside the grid");

    else if(X[0] > x_0+lx || X[1] > y_0+ly)
        throw runtime_error("Material point outside the grid");
        
    else
    {
        double dx_inv, dy_inv, a_0, b_0, a_1, b_1, da_0, db_0, da_1, db_1;
        array<array<double, 2>, 4> X_I;

        dx_inv = 1/dx; dy_inv = 1/dy;
        X_I = Cell_X(X);

        a_0 = (X_I[2][0] - X[0]) * dx_inv;
        a_1 = (X[0] - X_I[0][0]) * dx_inv;
        b_0 = (X_I[2][1] - X[1]) * dy_inv;
        b_1 = (X[1] - X_I[0][1]) * dy_inv;
            
        // da_0 = -X[0] * dx_inv;
        // da_1 =  X[0] * dx_inv;
        // db_0 = -X[1] * dy_inv;
        // db_1 =  X[1] * dy_inv;
        
        da_0 = -dx_inv;
        da_1 =  dx_inv;
        db_0 = -dy_inv;
        db_1 =  dy_inv;

        SHPC_mp[0] = a_0 * b_0; SHPC_mp[1]  = da_0 * b_0; SHPC_mp[2]  = a_0 * db_0;
        SHPC_mp[3] = a_1 * b_0; SHPC_mp[4]  = da_1 * b_0; SHPC_mp[5]  = a_1 * db_0;
        SHPC_mp[6] = a_1 * b_1; SHPC_mp[7]  = da_1 * b_1; SHPC_mp[8]  = a_1 * db_1;
        SHPC_mp[9] = a_0 * b_1; SHPC_mp[10] = da_0 * b_1; SHPC_mp[11] = a_0 * db_1;

        return SHPC_mp;
    }     
}

void MPM_Grid::ResetNC()
/*
    Reset the grids and set the values to zero

    Syntax: ResetNC()

    NC -> Nodal container, contains all computed nodal values

    The Nodal container has 11 columns with respective nodal values stored

    Nodal masses        | Nodal momentum (X)    | Nodal momentum (Y)        |

    Nodal forces (X)    | Nodal forces (Y)      | Nodal acceleration (X)    | Nodal acceleration (Y)|

    Nodal velocity (X)  | Nodal velocity (Y)    | Nodal displacecment (X)   | Nodal displacement (Y)|
*/
{
    for (size_t i = 0; i < (NoNodes) ; i++)
    {
        for (size_t j = 0; j < 11; j++)
        {
            NC[i][j] = 0;
        }        
    }
}

void MPM_Grid::SolveGrid(double dt, double mass_cutoff)
/*
        Solves for nodal acceleration, nodal velocity and nodal displacement and assembles to Nodal container

        Syntax - Solve_grid(NC,dt,mass_cutoff = 1e-10)

        NC - Nodal container of the Body

        dt - time step

        mass_cutoff - threshold mass value to avoid mathematical instability (optional argument)
*/
{
    NC.resize(NoNodes);
    double m_inv;
    double alpha = 0;      //numerical damping

    for (size_t i = 0; i < (NoNodes) ; i++)
    {
        if (NC[i][0] >= mass_cutoff)
        {
            m_inv = 1/NC[i][0];
            NC[i][5] = (NC[i][3] - (alpha * NC[i][1])) * m_inv;     //Nodal acceleration(X)
            NC[i][6] = (NC[i][4] - (alpha * NC[i][2])) * m_inv;     //Nodal acceleration(Y)

            NC[i][7] = (NC[i][1] * m_inv) + NC[i][5] * dt;          //Nodal velocity(X)
            NC[i][8] = (NC[i][2] * m_inv) + NC[i][6] * dt;          //Nodal velocity(Y)

            NC[i][9] = NC[i][7] * dt;                               //Nodal displacement(X)
            NC[i][10] = NC[i][8] * dt;                              //Nodal displacement(Y)
        }
    }
}

vector<array<double, 2>> MPM_Grid::Grid_X()
{
    vector<array<double, 2>> GX_I;      //Grid Cooedinates
    GX_I.resize(NoNodes);
    
    //Initialization
    for (size_t i = 0; i < (NoNodes); i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            GX_I[i][j] = 0;
        }
    }
        
    //Computation
    for (size_t i = 0; i < (nx+1); i++)
    {
        for (size_t j = 0; j < (ny+1); j++)
        {
            GX_I[(j*(nx+1))+i][0] = x_0 + (i*dx);
            GX_I[(j*(nx+1))+i][1] = y_0 + (j*dy);
        }
    }

    return GX_I;
}

void MPM_Grid::VTKGrid(string output_file_name, size_t step)
/*
    Prepares the vtk files for post processing

    Syntax VTKGrid("output_file_name", step)

    step - The iteration
*/
{
    //Declaring variables
    vector<array<double, 2>> GX_I;      //Grid Cooedinates 
    vector<array<size_t, 4>> GN_I;      //Grid connection
    array<double,2> X;                  //temporary material points 
    array<size_t,4> N_I;                //Local index
    
    GX_I.resize(NoNodes);
    GN_I.resize(nx*ny);

    //Initialization
    for (size_t i = 0; i < (NoNodes); i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            GX_I[i][j] = 0;
        }
    }

    for (size_t i = 0; i < (nx*ny); i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            GN_I[i][j] = 0;
        }
    }
    
    //Computation
    for (size_t i = 0; i < (nx+1); i++)
    {
        for (size_t j = 0; j < (ny+1); j++)
        {
            GX_I[(j*(nx+1))+i][0] = x_0 + (i*dx);
            GX_I[(j*(nx+1))+i][1] = y_0 + (j*dy);
        }
    }

    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            X = {(x_0 + (dx*0.5) + (i*dx)),(y_0 + (dy*0.5) + (j*dy))};
            N_I = Cell_I(X);
            GN_I[(i*ny)+j][0] = N_I[0]; GN_I[(i*ny)+j][1] = N_I[1]; GN_I[(i*ny)+j][2] = N_I[2]; GN_I[(i*ny)+j][3] = N_I[3];
        }
    }

    //Post processing Variables
    vector<string> ScalarPointDataNames;
    vector<vector<double> > ScalarPointData;
    vector<string> VectorPointDataNames;
    vector<vector<array<double, 2>> > VectorPointData;

    vector<double> M_I;
    vector<array<double, 2>> mv_I, f_I;
    M_I.resize(NoNodes); mv_I.resize(NoNodes); f_I.resize(NoNodes);

    //Post processing
    for (size_t j = 0; j < (NoNodes); j++)
    {
        M_I[j] = NC[j][0];          //Nodal masses
        mv_I[j][0] = NC[j][1];      //Nodal momentum
        mv_I[j][1] = NC[j][2];      //Nodal momentum
        f_I[j][0] = NC[j][3];       //Nodal forces
        f_I[j][1] = NC[j][4];       //Nodal forces
    }

    //Postprocessing data
    ScalarPointDataNames.push_back("Nodal_masses");
    ScalarPointData.push_back(M_I);
    //ScalarPointData = {M_I};
    VectorPointDataNames.push_back("Nodal_momentum");
    VectorPointData.push_back(mv_I);
    VectorPointDataNames.push_back("Nodal_forces");
    VectorPointData.push_back(f_I);
    
    vtk_export(output_file_name, step, GX_I, GN_I, ScalarPointDataNames, ScalarPointData, VectorPointDataNames, VectorPointData);

}