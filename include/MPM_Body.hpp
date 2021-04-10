#pragma once 

#include<iostream>
#include<array>
#include<vector>
#include<string>
#include<cmath>
#include<output_util.hpp>
#include<MPM_Grid.hpp>

using namespace std;

class MPM_Body
{    
    public:
    vector<array<double,19>> MPC;       //Material point container
    vector<array<size_t, 4>> IDC;       //Declaring Node container
    vector<array<double, 12>> SHPC;     //Declaring Shape function container
    size_t NoMP;                        //Number of material points

    public:
    //constructor
    MPM_Body(vector<array<double,3>> &MPV)
    /*    
        Creates and initialises the Material point container(MPC), Local nodal index container(IDC) and Shape function container(SHPC)

        MPV - An array containing the X and Y and volume per particle of every material points
        
        --------------------------------------------------------------

        The MPC has 19 columns with respective values of the material point stored

        X Coordinate                  | Y Coordinate                  |
    
        Velocity in X direction       | Velocity in Y direction       | Volume per particle           |
    
        Cauchy Stress along XX        | Cauchy Stress along XY        | Cauchy Stress along YY        | Cauchy Stress along ZZ        |

        Deformation Gradient along XX | Deformation Gradient along XY | Deformation Gradient along YX | Deformation Gradient along YY |

        Velocity Gradient along XX    | Velocity Gradient along XY    | Velocity Gradient along YX    | Velocity Gradient along YY    |

        Gravity along X               |  Gravity along Y              |

        ---------------------------------------------------------------
        
        IDC - Local nodal index container for material points 

        Node 1 | Node 2 | Node 3 | Node 4 |

        ----------------------------------------------------------------

        SHPC - Shape function container for material points

        Node 1 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |,

        Node 2 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |,

        Node 3 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |,

        Node 4 ->   Shape function | Derived shape function wrt X | Derived shape function wrt Y |

    */  
    {   
        MPV.shrink_to_fit();
        NoMP = MPV.size();
        MPC.resize(NoMP);   SHPC.resize(NoMP);  IDC.resize(NoMP);

        //Initialising the containers
        //Nodal Container
        for (size_t i = 0; i <  NoMP; i++)
        {
            for (size_t j = 0; j < 4; j++)
            {
                IDC[i][j] = 0.0;
            }
        }
        //Shape functoin Container
        for (size_t i = 0; i < NoMP; i++)
        {
            for (size_t j = 0; j < 12; j++)
            {
                SHPC[i][j] = 0.0;
            }
        }
        //Material point container
        for (size_t i = 0; i < NoMP; i++)
        {
            for (size_t j = 0; j < 19; j++)
            {
                MPC[i][j] = 0.0;
            }
        }

        //Loading data into Material container
        for (size_t i = 0; i < NoMP; i++)
        {
            MPC[i][0] = MPV[i][0];          //X Coordinates
            MPC[i][1] = MPV[i][1];          //Y Coordinates
            MPC[i][4] = MPV[i][2];          //Volume of each material point
            MPC[i][9] = 1; MPC[i][12] = 1;  //Deformation gradient
        }      
    }

    public:
    void setvelocity(vector<double> &vx, vector<double> &vy)
    /*
        Updates the velocity of the material points in the MPC

        Syntax: setvelocity(vx,vy)

        vx, vy - Velocity of particle in X and Y direction

    */
    {    
        vx.resize(NoMP);    vy.resize(NoMP);
        for (size_t i = 0; i < NoMP; i++)
        {
            MPC[i][2] = vx[i]; MPC[i][3] = vy[i];
        }
    }    

    public:
    array<array<double, 2>, 2> getdeformationgradient(size_t p)
    /*
        Returns the Deformation gradient of the pth material point

        Syntax: getDeformationGradient(p)

        p - Material point number (NOTE: Material point numbering starts from 0)

    */
    {   
        array<array<double, 2>, 2> F;                       //Declaration of deformation gradient
        F[0][0] = 0; F[1][1]=0;  F[1][0] = 0; F[0][1] = 0;  //Initialisation of deformation gradient

        F[0][0] = MPC[p][9]; F[0][1] = MPC[p][10]; F[1][0] = MPC[p][11]; F[1][1] = MPC[p][12];      //Retrieving deformation gradient of pth particle
        return F;
    }

    public:
    void ComputeMapping(MPM_Grid &Grid)
    /*
        Maps the material point parameters to the grid using shape functions and updates IDC and SHPC.

        Syntax: Compute_mapping(Grid)

        Grid - Computational Background Grid

    */
    {   
        //Declaration
        vector<array<double,2>> X;      //Coordinates of material points
        array<double,2> x;              //Coordinates of ith material point
        array<size_t,4> N_I;            //Local nodal index of ith material point
        array<double,12> SHP;           //Shape functions of ith material point

        X.resize(NoMP); IDC.resize(NoMP); SHPC.resize(NoMP);
        for (size_t i = 0; i < NoMP; i++)
        {
            X[i][0] = MPC[i][0];    X[i][1] = MPC[i][1];
        }
    
        for (size_t i = 0; i < NoMP; i++)
        {
            x[0] = X[i][0];     x[1] = X[i][1];
            N_I = Grid.Cell_I(x);
            
            IDC[i][0] = N_I[0];     IDC[i][1] = N_I[1];     IDC[i][2] = N_I[2];     IDC[i][3] = N_I[3];
        }
    
        for (size_t i = 0; i < NoMP; i++)
        {
            x[0] = X[i][0];     x[1] = X[i][1];
            SHP = Grid.SHP(x);

            for (size_t j = 0; j < 12; j++)
            {
                SHPC[i][j] = SHP[j];
            } 
        }
    }

    public:
    virtual void ComputeRHS(MPM_Grid &Grid) = 0;

    public:
    virtual array<array<double, 3>, 3> getstress(size_t i, MPM_Grid &Grid) = 0;

    public:
    virtual void Update(MPM_Grid &Grid, double dt) = 0;

    public:
    vector<double> KE(double rho)
    /*
        Returns the total Kinetic energy of the body

        Syntax: KE(rho)

        rho - Density of the body

    */
    {
        vector<double> KE;
        KE.resize(NoMP);
        double v;
        for (size_t i = 0; i < NoMP; i++)
        {
           v = sqrt ( (pow(MPC[i][2],2)) + (pow(MPC[i][3],2)) );
           KE[i] = 0.5 * v*v * rho*MPC[i][4];
        }      
        return KE;
    }
    
    public:
    vector<double> SE()
    /*
        Returns the total Strain energy of the body

        Syntax: SE()

    */
    {
        vector<double> SE;
        SE.resize(NoMP);
        array<array<double, 2>, 2> Sig, F, eps;
        double mul;
        
        for (size_t i = 0; i < 2; i++)
        {
        for (size_t j = 0; j < 2; j++)
            {
                Sig[i][j] = 0; F[i][j] = 0; eps[i][j] = 0;           
            }        
        }

        double v;
        for (size_t i = 0; i < NoMP; i++)
        {
           mul = 0;
           F = getdeformationgradient(i);

           eps[0][0] = F[0][0]-1; eps[1][1] = F[1][1]-1;
           eps[0][1] = F[0][1]; eps[1][0] = F[1][0];

            Sig[0][0] = MPC[i][5]; Sig[1][1] = MPC[i][7]; Sig[2][2] = MPC[i][8];
            Sig[0][1] = MPC[i][6]; Sig[1][0] = MPC[i][6];

           for (size_t l = 0; l < 2; l++)
           {
               for (size_t m = 0; m < 2; m++)
               {
                   mul = mul + (Sig[l][m] * eps[l][m]); 
               }               
           }

           SE[i] = 0.5 * mul * MPC[i][4];         
        }

        return SE;
    }
    public:
    void VTKMaterialpoints(string output_file_name, size_t step)
    /*
        Prepares the vtk files for post processing

        Syntax VTKMaterialpoints("output_file_name", step)

        step - The iteration
    */
    {
        vector<array<double, 2>> X_mp;      //Material point coordinates
        X_mp.resize(NoMP);
        //Initialization
        for (size_t i = 0; i < NoMP; i++)
        {
            for (size_t j = 0; j < 2; j++)
            {
                X_mp[i][j] = 0;
            }
        }

        //Computation
        for (size_t i = 0; i < NoMP; i++)
        {
            X_mp[i][0] = MPC[i][0];
            X_mp[i][1] = MPC[i][1];
        }

        //Declaring post processing Variables
        vector<string> ScalarPointDataNames;
        vector<vector<double> > ScalarPointData;
        vector<string> VectorPointDataNames;
        vector<vector<array<double, 2>> > VectorPointData;

        vector<double> Vol, VonMisses;      //Volume of material points
        vector<array<double, 2>> v;         //velocity of material points
        Vol.resize(NoMP); VonMisses.resize(NoMP); v.resize(NoMP);

        //Post processing
        for (size_t j = 0; j < (NoMP); j++)
        {
            VonMisses[j] = sqrt ( (pow(MPC[j][5] - MPC[j][7],2)) + (pow(MPC[j][7] - MPC[j][8],2)) +
                                  (pow(MPC[j][5] - MPC[j][8],2)) + ((6 * pow(MPC[j][6],2))/2) );

            Vol[j] = MPC[j][4];
            v[j][0] = MPC[j][2];
            v[j][1] = MPC[j][3];
        }
        ScalarPointDataNames.push_back("Volume");
        ScalarPointData.push_back(Vol);
        ScalarPointDataNames.push_back("VonMisses");
        ScalarPointData.push_back(VonMisses);
        VectorPointDataNames.push_back("Velocity");
        VectorPointData.push_back(v);

        //Exporting ot .vtk files
        vtk_export(output_file_name, step, X_mp, ScalarPointDataNames, ScalarPointData, VectorPointDataNames, VectorPointData);
    }
    
    public:
    //destructor
    ~MPM_Body(){}

};