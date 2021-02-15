"""
Bench-marking Two Disks

Problem description: 
Young's modulus = 1e4, Density = 1000, Poisson ratio = 0.3
No. of Material points = 10
No. of cells  = 10
Gravity = 9.81

"""
import os,sys,shutil
import numpy as np
import matplotlib.pyplot as plt
import time
from mpm_circ import mpm_circ
sys.path.append(os.getcwd()+'/build')
import MPM_2D
from tqdm import tqdm, trange

#Defining the Disks
Emod = 1e4; Density = 1000; Poisson = 0.3               #Material parameters
MPV_1 = mpm_circ([0.2,0.2], 0.2, [10,10], 2)
MPV_2 = mpm_circ([0.8,0.8], 0.2, [10,10], 2)

Disk_1 = MPM_2D.MPM_Solid(Emod, Density, Poisson, MPV_1)
Disk_1.MPC()[:,2] = 0.1; Disk_1.MPC()[:,3] = 0.1

Disk_2 = MPM_2D.MPM_Solid(Emod, Density, Poisson, MPV_2)
Disk_2.MPC()[:,2] = -0.1; Disk_2.MPC()[:,3] = -0.1 

#Defining the grid
x_0 = 0; y_0 = 0;                                       #Origin of X and Y Axis
lx = 1.0; ly = 1.0;                                     #Length of grid in X and Y direction
nx = 20; ny = 20;                                       #Number of cells in X and Y direction
Grid = MPM_2D.MPM_Grid(L=[lx, ly], N=[nx, ny], O=[x_0, y_0])

#Prepare postprocessing
target_directory = os.path.join(os.getcwd(), "TwoDisks_Benchmarking")

#Create a new directory to store the plots (delete when there)
if os.path.exists(target_directory) and os.path.isdir(target_directory):
    shutil.rmtree(target_directory)
os.mkdir(target_directory)
#Grid.VTKGrid("/home/karthik/cpp/mpm_2d_cpp/TwoDisks_Benchmarking/DisksGrid", 0)

TwodisksSolver = MPM_2D.MPM_Solver()                    #Defining solver
Bodies = [Disk_1, Disk_2]

#Running through time steps
dt = 1e-4                                               #Time interval
T  = 3.0                                                #Total time period
NoS = int(T/dt)
noMicroSteps = int(100)                                  #Incremental time step for solver
noMacroSteps = int(NoS/noMicroSteps)                    #Number of vtk files
start_time = time.time()
for step in trange (noMacroSteps):
    TwodisksSolver.Solve(Grid, Bodies, noMicroSteps, dt)
    Grid.VTKGrid("/home/karthik/MPM_2D/TwoDisks_Benchmarking/DisksGrid", step)
    Disk_1.VTKMaterialpoints("/home/karthik/MPM_2D/TwoDisks_Benchmarking/Disk_1", step)
    Disk_2.VTKMaterialpoints("/home/karthik/MPM_2D/TwoDisks_Benchmarking/Disk_2", step)
print("--- %s seconds ---" % (time.time() - start_time))