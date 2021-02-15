#include<MPM_Body.hpp>

class MPM_Fluid : public MPM_Body
{   
    public:
    double mu, rho, kappa, gamma;

    public:
    //constructor
    MPM_Fluid(double visc, double density, double bulk_mod, vector<array<double,3>> &MPV, double isen_exp = 7) : MPM_Body(MPV)
    /*    
        Initialises the Fluid Body

        Syntax: MPM_Fluid(visc, density, bulk_mod)

        visc - Viscosity of the fluid

        rho - Density of the fluid

        bulk_mod = Bulk Modulus of the fluid

        gammma = isentropic exponent

    */ 
    {
        mu = visc;      
        rho = density;
        kappa = bulk_mod;
        gamma = isen_exp;
    
    }   
    
    public:
    void ComputeRHS(MPM_Grid &Grid);

    public:
    array<array<double, 3>, 3> getstress(size_t i, MPM_Grid &Grid);

    public:
    void Update(MPM_Grid &Grid, double dt);

    public:
    //destructor
    ~MPM_Fluid(){}

};