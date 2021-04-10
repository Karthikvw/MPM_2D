"""
Bench-marking 1D Spring

Problem description: 
Young's modulus = 1e5, Density = 1500, Poisson ratio = 0
No. of Material points = 10
No. of cells  = 10
Gravity = 9.81

"""
import os,sys,shutil
import numpy as np
import time
import matplotlib.pyplot as plt
from mpm_rect import mpm_rect
sys.path.append(os.getcwd()+'/build')
import MPM_2D
from tqdm import tqdm, trange

#Defining the Spring
Emod = 1e5; Density = 1500; Poisson = 0.3                                   #Material parameters
MPV = mpm_rect([0.06,0.12], [0.06,0.06], [9,1000],1  )
Truss = MPM_2D.MPM_NHolid(Emod, Density, Poisson, MPV)
Truss.MPC()[:,18] = -1400                                                  #Body force(gravity)

#Defining the grid
x_0 = 0; y_0 = 0;                                                           #Origin of X and Y Axis
lx = 0.18; ly = 0.18;                                                     #Length of grid in X and Y direction
nx = 3; ny = 3;                                                             #Number of cells in X and Y direction
Grid = MPM_2D.MPM_Grid(L=[lx, ly], N=[nx, ny], O=[x_0, y_0])

#Prepare postprocessing
target_directory = os.path.join(os.getcwd(), "Truss_Benchmarking")

#Create a new directory to store the plots (delete when there)
if os.path.exists(target_directory) and os.path.isdir(target_directory):
    shutil.rmtree(target_directory)
os.mkdir(target_directory)
#Grid.VTKGrid("/home/karthik/cpp/mpm_2d_cpp/Truss_Benchmarking/TrussGrid", 0)

TrussSolver = MPM_2D.MPM_Solver()                                          #Defining solver
alpha = 25
TrussSolver.addBC(lambda X: True if (X[1] >= 0.18) else False, vy=0.0)     #Boundary conditions

#Running through time steps
dt = 1e-4                                                                   #Time interval
T  = 2.0                                                                    #Total time period
NoS = int(T/dt)
noMicroSteps = int(10)                                                     #Incremental time step for solver
noMacroSteps = int(NoS/noMicroSteps)                                        #Number of vtk files
start_time = time.time()
for step in trange (noMacroSteps):
    TrussSolver.Solve(Grid, [Truss], noMicroSteps, dt, alpha)
    Grid.VTKGrid("/home/karthik/MPM_2D/Truss_Benchmarking/TrussGrid", step)
    Truss.VTKMaterialpoints("/home/karthik/MPM_2D/Truss_Benchmarking/Truss", step)
print("--- %s seconds ---" % (time.time() - start_time))