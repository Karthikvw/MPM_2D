#include<MPM_Fluid.hpp>
#include<cmath>

void MPM_Fluid::ComputeRHS(MPM_Grid &Grid)
/*
    Computes nodal mass vector, nodal momentum vector and nodal force vector and assembles it to Nodal container

        Syntax: Compute_RHS(Grid)

        Grid -> Computautional Background Grid

*/
{
    //Declaring Variables
    vector<array<double,4>> SHP, dSHPx, dSHPy, Sig;      //Shape functions, derived shape functions and Cauchy stresses 
    vector<double> M, mvx, mvy;                          //Mass and momentum vector
    array<array<double, 2>, 2> F;                        //Particle deformation gradient
    double J;                                            //Jacobian
    
    //Resizing all declared variables
    SHP.resize(NoMP); dSHPx.resize(NoMP); dSHPy.resize(NoMP);   
    M.resize(NoMP); mvx.resize(NoMP); mvy.resize(NoMP);
    Sig.resize(NoMP);

    //Initialising variables
    for (size_t i = 0; i < NoMP; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            SHP[i][j] = 0; dSHPx[i][j]=0; dSHPy[i][j] = 0; Sig[i][j] = 0;
        }        
    }

    for (size_t i = 0; i < NoMP; i++)
    {
        M[i] = 0; mvx[i] = 0; mvy[i] = 0;
    }    

    // //Retrieving the data from the containers
    for (size_t i = 0; i < NoMP; i++)
    {
        SHP[i][0] = SHPC[i][0]; SHP[i][1] = SHPC[i][3]; SHP[i][2] = SHPC[i][6]; SHP[i][3] = SHPC[i][9];             //Shape functions
        dSHPx[i][0] = SHPC[i][1]; dSHPx[i][1] = SHPC[i][4]; dSHPx[i][2] = SHPC[i][7]; dSHPx[i][3] = SHPC[i][10];    //Derived shape functions(X)
        dSHPy[i][0] = SHPC[i][2]; dSHPy[i][1] = SHPC[i][5]; dSHPy[i][2] = SHPC[i][8]; dSHPy[i][3] = SHPC[i][11];    //Derived shape functions(X)
        M[i] = rho * MPC[i][4];                                                                                     //Mass Vector
        mvx[i] = rho * MPC[i][4] * MPC[i][2];                                                                       //Momentum Vector(X)
        mvy[i] = rho * MPC[i][4] * MPC[i][3];                                                                       //Momentum Vector(X)
        Sig[i][0] = MPC[i][5]; Sig[i][1] = MPC[i][6]; Sig[i][2] = MPC[i][6]; Sig[i][3] = MPC[i][7];                 //Cauchy Stresses
    }

    //Corresponding variables for ith material point
    array<double,4> m_I_mp, mvx_I_mp, mvy_I_mp, fx_I_mp, fy_I_mp;
    array<double,4> SHP_mp, dSHPx_mp, dSHPy_mp, Sig_mp;
    array<size_t,4> ID_mp;
    for (size_t i = 0; i < NoMP; i++)      //Looping through all material points
    {
        for (size_t j = 0; j < 4; j++)     //Looping through 4 nodes of ith material point
        {
            //Initialising parameters
            m_I_mp[j]  = 0;     SHP_mp[j] = 0;
            mvx_I_mp[j] = 0;    dSHPx_mp[j] = 0;
            mvy_I_mp[j] = 0;    dSHPy_mp[j] = 0;
            fx_I_mp[j] = 0;     ID_mp[j]  = 0;
            fy_I_mp[j] = 0;     Sig_mp[j] = 0;

            //Retrieving for nth material point
            SHP_mp[j]   = SHP[i][j];
            dSHPx_mp[j] = dSHPx[i][j];
            dSHPy_mp[j] = dSHPy[i][j];
            ID_mp[j]    = IDC[i][j];
            Sig_mp[j]   = Sig[i][j];
        }

        for (size_t j = 0; j < 4; j++)     //Looping through 4 nodes of each material point
        {
            //Computing and assembling Nodal Masses
            m_I_mp[j] = M[i] * SHP_mp[j];
            Grid.NC[ID_mp[j]][0] +=  m_I_mp[j]; 

            //Computing and assembling nodal momentum vector
            mvx_I_mp[j] = mvx[i] * SHP_mp[j];
            mvy_I_mp[j] = mvy[i] * SHP_mp[j];
            Grid.NC[ID_mp[j]][1] += mvx_I_mp[j];
            Grid.NC[ID_mp[j]][2] += mvy_I_mp[j];

            //Computing and assembling nodal force vector
            F = getdeformationgradient(i);
            J = (F[0][0] * F[1][1]) - (F[1][0] * F[0][1]);

            fx_I_mp[j] = MPC[i][17]*M[i]*SHP_mp[j] - dSHPx_mp[j]*Sig_mp[0]*J*MPC[i][4] - dSHPy_mp[j]*Sig_mp[1]*J*MPC[i][4];
            fy_I_mp[j] = MPC[i][18]*M[i]*SHP_mp[j] - dSHPx_mp[j]*Sig_mp[1]*J*MPC[i][4] - dSHPy_mp[j]*Sig_mp[3]*J*MPC[i][4];
            Grid.NC[ID_mp[j]][3] += fx_I_mp[j];
            Grid.NC[ID_mp[j]][4] += fy_I_mp[j];
        }    
    }    
}

array<array<double, 3>, 3> MPM_Fluid::getstress(size_t i, MPM_Grid &Grid)
/*
    Get the Cauchy stresses for the ith material point

    Syntax: getstress(i, Grid)

    i -> material point index

    Grid -> Computational Background Grid

*/

{
    //Declaring Variables
    array<array<double, 3>, 3> Sig;                         //Cauchy stresses
    array<array<double, 3> ,3> D, devD, I;                  //Strain rate, Deviotoric part of D, Identity tensor
    array<array<double, 2>, 2> f, L;                        //Deformation Gradient, Gradient of the velocity    
    double trD, p;                                          //Trace of D, Pressure term
    double J_inv, J;                                        //Jacobian and inverse of Jacobian
    double theta;                                           //Volumetric dilation

    array<double ,4> ax_I, ay_I, vx_I, vy_I;                //Nodal acceleration and nodal velocity of ith material point
    array<size_t,4> ID_mp;                                  //Local nodal coordinates of ith materialpoint

    //Initialising the tensors
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            Sig[i][j] = 0; D[i][j] = 0; devD[i][j] = 0; I[i][j] = 0;
        }        
    }
    I[0][0] = 1; I[1][1] = 1; I[2][2] = 1;

    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            f[i][j] = 0; L[i][j] = 0;
        }        
    }    

    for (size_t j = 0; j < 4; j++)      //Looping through 4 nodes of ith material point
    {
        ID_mp[j]    = IDC[i][j];
        ax_I[j] = Grid.NC[ID_mp[j]][5];                      //Nodal acceleration(X) of ith material point
        ay_I[j] = Grid.NC[ID_mp[j]][6];                      //Nodal acceleration(Y) of ith material point
        vx_I[j] = Grid.NC[ID_mp[j]][7];                      //Nodal velocity(X) of ith material point
        vy_I[j] = Grid.NC[ID_mp[j]][8];                      //Nodal velocity(Y) of ith material point
    }

     //Computing gradient of the velocity of the material point
     L[0][0] = SHPC[i][1]*vx_I[0] + SHPC[i][4]*vx_I[1] + SHPC[i][7]*vx_I[2] + SHPC[i][10]*vx_I[3];
     L[0][1] = SHPC[i][2]*vx_I[0] + SHPC[i][5]*vx_I[1] + SHPC[i][8]*vx_I[2] + SHPC[i][11]*vx_I[3];
     L[1][0] = SHPC[i][1]*vy_I[0] + SHPC[i][4]*vy_I[1] + SHPC[i][7]*vy_I[2] + SHPC[i][10]*vy_I[3];
     L[1][1] = SHPC[i][2]*vy_I[0] + SHPC[i][5]*vy_I[1] + SHPC[i][8]*vy_I[2] + SHPC[i][11]*vy_I[3];

     D[0][0] = 0.5*(L[0][0] + L[0][0]); D[0][1] = 0.5*(L[0][1] + L[1][0]); 
     D[1][0] = 0.5*(L[1][0] + L[0][1]); D[1][1] = 0.5*(L[1][1] + L[1][1]);
     trD = D[0][0] + D[1][1] + D[2][2];
     theta = trD/3;

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            devD[i][j] = D[i][j] - ((trD/3) * I[i][j]);
        }        
    }

     f = getdeformationgradient(i);
     J = (f[0][0] * f[1][1]) - (f[1][0] * f[0][1]);
     J_inv = 1/((f[0][0] * f[1][1]) - (f[1][0] * f[0][1]));
     p = kappa * ((pow(J_inv, gamma)) - 1);

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            Sig[i][j] = (2*mu*devD[i][j]) - (p*I[i][j]);
        }        
    }
    
    return Sig;
}

void MPM_Fluid::Update(MPM_Grid &Grid, double dt)
/*
       Updates the material point parameters using USL Explicit method

        Syntax: Update(Grid,dt)

        Grid - Computational Background Grid

        dt - time increment
*/
{
     //Declaring Variables
     vector<array<double,4>> SHP, dSHPx, dSHPy;           //Shape functions, derived shape functions and Cauchy stresses
     array<double ,4> ax_I, ay_I, vx_I, vy_I;             //Nodal acceleration and nodal velocity of ith material point
     array<array<double, 3>, 3> Sig;                      //Cauchy stresses    
    
     array<double,4> SHP_mp, dSHPx_mp, dSHPy_mp;          //Shape functions, derived shape functions of ith material point
     array<size_t,4> ID_mp;                                  //Local nodal coordinates of ith materialpoint
     double ax_mp, ay_mp, vhx_mp, vhy_mp;                 //Acceleration and velocity of ith material point
     array<array<double, 2> ,2> F_mp, f_mp, gvh_mp, delF; //Deformation gradient, velocity gradient, Transformation matrix of ith material point
     array<array<size_t, 2> ,2> I;                           //Identity matrix
    
     //Resizing all declared variables
     IDC.resize(NoMP); SHPC.resize(NoMP);
     SHP.resize(NoMP); dSHPx.resize(NoMP); dSHPy.resize(NoMP);

     //Initialising variables
     for (size_t i = 0; i < NoMP; i++)
     {
         for (size_t j = 0; j < 4; j++)
         {
             SHP[i][j] = 0; dSHPx[i][j]=0; dSHPy[i][j] = 0;
         }        
     }

     for (size_t i = 0; i < 3; i++)
         {
         for (size_t j = 0; j < 3; j++)
             {
                 Sig[i][j] = 0;            
             }        
         }

     I[0][0] = 1; I[0][1] = 0; I[1][0] = 0; I[1][1] = 1;

     //Retrieving the data from the containers
     for (size_t i = 0; i < NoMP; i++)
     {
         SHP[i][0] = SHPC[i][0]; SHP[i][1] = SHPC[i][3]; SHP[i][2] = SHPC[i][6]; SHP[i][3] = SHPC[i][9];             //Shape functions
         dSHPx[i][0] = SHPC[i][1]; dSHPx[i][1] = SHPC[i][4]; dSHPx[i][2] = SHPC[i][7]; dSHPx[i][3] = SHPC[i][10];    //Derived shape functions(X)
         dSHPy[i][0] = SHPC[i][2]; dSHPy[i][1] = SHPC[i][5]; dSHPy[i][2] = SHPC[i][8]; dSHPy[i][3] = SHPC[i][11];    //Derived shape functions(X)
     }

     for (size_t i = 0; i < NoMP; i++)      //Looping through all material points
     {
        for (size_t p = 0; p < 2; p++)
        {
            for (size_t q = 0; q < 2; q++)
            {
                 F_mp[p][q] = 0; f_mp[p][q] = 0; gvh_mp[p][q]=0; delF[p][q] = 0;
            }        
        }

         for (size_t j = 0; j < 4; j++)     //Looping through 4 nodes of ith material point
         {
             //Initialising parameters
             SHP_mp[j] = 0;          ax_I[j] = 0;
             dSHPx_mp[j] = 0;        ay_I[j] = 0;
             dSHPy_mp[j] = 0;        vx_I[j] = 0;
             ID_mp[j]  = 0;          vy_I[j] = 0;

             //Retrieving for nth material point
             SHP_mp[j]   = SHP[i][j];
             dSHPx_mp[j] = dSHPx[i][j];
             dSHPy_mp[j] = dSHPy[i][j];
             ID_mp[j]    = IDC[i][j];
             ax_I[j] = Grid.NC[ID_mp[j]][5];                      //Nodal acceleration(X) of ith material point
             ay_I[j] = Grid.NC[ID_mp[j]][6];                      //Nodal acceleration(Y) of ith material point
             vx_I[j] = Grid.NC[ID_mp[j]][7];                      //Nodal velocity(X) of ith material point
             vy_I[j] = Grid.NC[ID_mp[j]][8];                      //Nodal velocity(Y) of ith material point
         }        

         //Computing acceleration of the material point along X and Y
         ax_mp = SHP_mp[0]*ax_I[0] + SHP_mp[1]*ax_I[1] + SHP_mp[2]*ax_I[2] + SHP_mp[3]*ax_I[3];
         ay_mp = SHP_mp[0]*ay_I[0] + SHP_mp[1]*ay_I[1] + SHP_mp[2]*ay_I[2] + SHP_mp[3]*ay_I[3];
            
         //Computing velocity of the material point along X and Y
         vhx_mp = SHP_mp[0]*vx_I[0] + SHP_mp[1]*vx_I[1] + SHP_mp[2]*vx_I[2] + SHP_mp[3]*vx_I[3];
         vhy_mp = SHP_mp[0]*vy_I[0] + SHP_mp[1]*vy_I[1] + SHP_mp[2]*vy_I[2] + SHP_mp[3]*vy_I[3];
            
         //Computing gradient of the velocity of the material point
         gvh_mp[0][0] = dSHPx_mp[0]*vx_I[0] + dSHPx_mp[1]*vx_I[1] + dSHPx_mp[2]*vx_I[2] + dSHPx_mp[3]*vx_I[3];
         gvh_mp[0][1] = dSHPy_mp[0]*vx_I[0] + dSHPy_mp[1]*vx_I[1] + dSHPy_mp[2]*vx_I[2] + dSHPy_mp[3]*vx_I[3];
         gvh_mp[1][0] = dSHPx_mp[0]*vy_I[0] + dSHPx_mp[1]*vy_I[1] + dSHPx_mp[2]*vy_I[2] + dSHPx_mp[3]*vy_I[3];
         gvh_mp[1][1] = dSHPy_mp[0]*vy_I[0] + dSHPy_mp[1]*vy_I[1] + dSHPy_mp[2]*vy_I[2] + dSHPy_mp[3]*vy_I[3];

         //Updating the position of the material point
         MPC[i][0] += vhx_mp * dt;
         MPC[i][1] += vhy_mp * dt;

         //Updating the velocity of the material point
         MPC[i][2] += ax_mp * dt;
         MPC[i][3] += ay_mp * dt;

         //Computing and Updating the deformation gradient of the material point
         F_mp = getdeformationgradient(i);
         for (size_t k = 0; k < 2; k++)
         {  
             for (size_t l = 0; l < 2; l++)
             {
                delF[k][l] = I[k][l] + (gvh_mp[k][l] * dt);
             }        
         }

         for (size_t m = 0; m < 2; m++)
         {
             for (size_t n = 0; n < 2; n++)
             {  
                 for (size_t o = 0; o < 2; o++)
                 {
                     f_mp[m][n] = f_mp[m][n] + (delF[m][o] * F_mp[o][n]);
                 }                 
             }        
         }
         MPC[i][9] = f_mp[0][0]; MPC[i][10] = f_mp[0][1]; MPC[i][11] = f_mp[1][0]; MPC[i][12] = f_mp[1][1];

         //Computing and updating the Cauchy stresses of the material point
         Sig = MPM_Fluid::getstress(i, Grid);
         MPC[i][5] = Sig[0][0]; MPC[i][6] = Sig[0][1]; MPC[i][7] = Sig[1][1]; MPC[i][8] = Sig[2][2];

     }

}