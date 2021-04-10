def mpm_rect(Origin, dim, Noc, Pattern):

    """
    The function discretises a rectangle and returns the material ponits with available two patterns only.

    Input Parameters:

    Origin: Bottom left coordinate of the rectangle

    dim: Dimensions of the rectangle in X and Y direction

    Noc: Number of cells in X and Y direction

    Pattern: Desired choice of dicretisation (Please enter only 1 or 2)

    Pattern 1: One material point at the centre of each cell

    Pattern 2: Four material points at each cell distributed based on Guass Legrandge quadrature rule

    Example: MP = mpm_rect([1,2], [6,4], [6,4],1)


    """

    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    #Conversion of data types and assigning variables
    Origin = np.array(Origin)
    x_0 = np.float(Origin[0]); y_0 = np.float(Origin[1]);
    lx = np.float(dim[0]); ly = np.float(dim[1]);
    nx = np.int(Noc[0]); ny = np.int(Noc[1]);
    
    dx = lx/nx; dy = ly/ny;             #Grid spacing in X and Y direction
    Vol = lx*ly                         #Total volume of the rectangle
   
    if (Pattern == 1):
        px = 0; py = 0;

        for i in range (nx):
            for j in range (ny):        
                px = np.append(px,(x_0+dx*i+x_0+dx*(i+1))/2)
                py = np.append(py,(y_0+dy*j+y_0+dy*(j+1))/2)
                
        px=px[1:(nx*ny)+1]; py=py[1:(nx*ny)+1]
        dVol = np.ones(nx*ny)*(Vol/(nx*ny))                                                        #Volume of each material point
        M=np.array((px,py,dVol))                                                                   #Array of Material points and dvol
        M=M.transpose() 
        return M

    if (Pattern == 2):
        w = 1/math.sqrt(3)                                                                         #Weights for Gaussian points
        px = 0; py=0;

        for i in range (nx):
            for j in range (ny):        
                px = np.append(px,(x_0+dx*i+x_0+dx*(i+1))/2 - dx*w*0.5)
                py = np.append(py,(y_0+dy*j+y_0+dy*(j+1))/2 - dy*w*0.5)
                
                px = np.append(px,(x_0+dx*i+x_0+dx*(i+1))/2 + dx*w*0.5)
                py = np.append(py,(y_0+dy*j+y_0+dy*(j+1))/2 + dy*w*0.5)
                
                px = np.append(px,(x_0+dx*i+x_0+dx*(i+1))/2 - dx*w*0.5)
                py = np.append(py,(y_0+dy*j+y_0+dy*(j+1))/2 + dy*w*0.5)
                
                px = np.append(px,(x_0+dx*i+x_0+dx*(i+1))/2 + dx*w*0.5)
                py = np.append(py,(y_0+dy*j+y_0+dy*(j+1))/2 - dy*w*0.5)
                
        px=px[1:(nx*ny*4)+1]; py=py[1:(nx*ny*4)+1] 
        dVol = np.ones(nx*ny*4)*(Vol/(len(px)))                                                    #Volume of each material point
        M=np.array((px,py,dVol))                                                                   #Array of Material points and dvol
        M=M.transpose()
        return M

    if (Pattern != 1 and Pattern != 2):
        print("Valid entries are only 1 and 2")
        print("*** Nothing returned ***")
