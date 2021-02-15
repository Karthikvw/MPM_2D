#include<MPM_Body.hpp>

class MPM_NHSolid : public MPM_Body
{
    public:
    double Emod, rho, nu;

    public:
    //constructor
    MPM_NHSolid(double E, double density, double pois_ratio, vector<array<double,3>> &MPV) : MPM_Body(MPV)
    /*    
        Initialises the Solid Body

        Syntax: MPM_Solid(E, rho, nu)

        E - Young's Modulus of the material

        rho - Density of the material

        nu - Poisson's ratio of the material

    */  
    {
        Emod = E;
        rho = density;
        nu = pois_ratio;
    }

    public:
    void ComputeRHS(MPM_Grid &Grid);

    public:
    array<array<double, 3>, 3> getstress(size_t i, MPM_Grid &Grid);

    public:
    void Update(MPM_Grid &Grid, double dt);

    public:
    //destructor
    ~MPM_NHSolid(){}
};