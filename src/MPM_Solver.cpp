#include<MPM_Solver.hpp>

void MPM_Solver::Solve(MPM_Grid* Grid, vector<MPM_Body*> Bodies, size_t NoS, double dt)
{
    double time = 0;                    //Time initilization
    for (size_t step=0; step < NoS; step++)
    {        
        Grid -> ResetNC();                                              //Reset Grid
        for (auto &Body : Bodies) Body -> ComputeMapping(*Grid);        //Mapping values to the grid
        for (auto &Body : Bodies) Body -> ComputeRHS(*Grid);            //Assembling values to the Nodal container
        Grid -> SolveGrid(dt);
        //Applying BC
        vector<array<double, 2>> XI = Grid -> Grid_X();                 //Get all nodal coordinates
        for (size_t n = 0; n < (Grid -> NoNodes); n++)
        {
            array<double, 2> X = XI[n];                                 //Get the nth nodal coordinate
            for (auto &BC : BoundaryConditions)                         //Loop over all possible boundary condition active at this node
            {
                if (BC.Criterion(X))
                {
                    if (BC.active_vx) {Grid -> NC[n][5] = 0.0; Grid -> NC[n][7] = BC.value_vx;}
                    if (BC.active_vy) {Grid -> NC[n][6] = 0.0; Grid -> NC[n][8] = BC.value_vy;}
                }
            }
        }

    // for (size_t i = 0; i < Bodies[0]-> NoMP ; i++)
    // {
    //     for (size_t j = 0; j < Bodies[0]-> MPC[0].size() ; j++)
    //     {
    //         cout<<Bodies[0]-> MPC[i][j]<<"\t";
    //     }
    //     cout<<endl;                
    // }
     for (auto &Body : Bodies) Body -> Update(*Grid, dt);               //Remapping values to particles

    // for (size_t i = 0; i < Bodies[0]-> NoMP ; i++)
    // {
    //     for (size_t j = 0; j < Bodies[0]-> MPC[0].size() ; j++)
    //     {
    //         cout<<Bodies[0]-> MPC[i][j]<<"\t";
    //     }
    //     cout<<endl;                
    // }
     
     time += dt;                                                        //Advance in time step

    }
}
