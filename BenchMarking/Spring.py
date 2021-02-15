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

#Defining the Spring
Emod = 1e5; Density = 1500; Poisson = 0.0                                   #Material parameters
MPV = mpm_rect([0,0.8], [0.02,0.2], [1,10],1  )
Spring = MPM_2D.MPM_Solid(Emod, Density, Poisson, MPV)
Spring.MPC()[:,18] = -9.81                                                  #Body force(gravity)

#Defining the grid
x_0 = 0; y_0 = 0;                                                           #Origin of X and Y Axis
lx = 0.02; ly = 1.0;                                                        #Length of grid in X and Y direction
nx = 1; ny = 21;                                                            #Number of cells in X and Y direction
Grid = MPM_2D.MPM_Grid(L=[lx, ly], N=[nx, ny], O=[x_0, y_0])

#Prepare postprocessing
target_directory = os.path.join(os.getcwd(), "Spring_Benchmarking")

#Create a new directory to store the plots (delete when there)
if os.path.exists(target_directory) and os.path.isdir(target_directory):
    shutil.rmtree(target_directory)
os.mkdir(target_directory)
#Grid.VTKGrid("/home/karthik/cpp/mpm_2d_cpp/Spring_Benchmarking/SpringGrid", 0)

#Running through time steps
t = 0                                                                       #Start time
T = 2.0; dt = 1e-4                                                          #Total time and time step
NoS = int(T/dt)                                                             #Number of time steps

SpringSolver = MPM_2D.MPM_Solver()                                          #Defining solver
SpringSolver.addBC(lambda X: True if (X[1] >= 1.0) else False, vy=0.0)      #Boundary conditions

#Running through time steps
dt = 1e-4                                                                   #Time interval
T  = 2.0                                                                    #Total time period
NoS = int(T/dt)
noMicroSteps = int(100)                                                     #Incremental time step for solver
noMacroSteps = int(NoS/noMicroSteps)                                        #Number of vtk files
start_time = time.time()
for step in range (noMacroSteps):
    SpringSolver.Solve(Grid, [Spring], noMicroSteps, dt)
    Grid.VTKGrid("/home/karthik/MPM_2D/Spring_Benchmarking/SpringGrid", step)
    Spring.VTKMaterialpoints("/home/karthik/MPM_2D/Spring_Benchmarking/spring_2d", step)
print("--- %s seconds ---" % (time.time() - start_time))