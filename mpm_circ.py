def mpm_circ(centre, rad, Noc, Pattern):

    """
    The function discretices a circle and returns the material points with available two patterns only.

    Input Parameters:

    Centre: Coordinate of the circle's centre (float)

    rad: Radius of the circle (float)

    Noc: Number of cells in X and Y direction (int)

    Pattern: Desired choice of dicretisation (Please enter only 1 or 2)
    
    (a material point is created only if the point lies inside the circle)

    Pattern 1: One material point at the centre of each cell

    Pattern 2: Four material points at each cell distributed based on Guass Legrandge quadrature rule

    Example: MP = mpm_circ([1,2], 3, [3,3],2)


    """

    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Circle
    
    centre = np.array(centre)
    x_c = np.float(centre[0]); y_c = np.float(centre[1]);
    r = np.float(rad);
    nx = np.int(Noc[0]); ny = np.int(Noc[1]);

    x_0 = x_c-r; y_0= y_c-r;                #Origin of enclosing square
    l = 2*r;                                #Length of the square
    dx = l/nx; dy = l/ny;                   #Grid spacing in X and Y direction

    Vol = math.pi*r*r;                      #Total Volume

    if (Pattern == 1):
        #fig, ax = plt.subplots(1,1,figsize = (5,5))
        #ax.set_xlim(x_0,l+x_0); ax.set_ylim(y_0,l+y_0);
        px = 0; py = 0

        for i in range (nx):
            #ax.plot([x_0+dx*i,x_0+dx*i],[y_0,l+y_0], c=(0.9, 0.9, 0.9, 1.0) )              #Plot X lines
            for j in range (ny):        
                #ax.plot([x_0,l+x_0],[y_0+dy*j,y_0+dy*j], c=(0.9, 0.9, 0.9, 1.0) )          #Plot Y lines
                a = (x_0+dx*i+x_0+dx*(i+1))/2; b = (y_0+dy*j+y_0+dy*(j+1))/2;
                if ( math.sqrt(((a-x_c)**2) + ((b-y_c)**2)) <= r   ):
                    px = np.append(px,a)
                    py = np.append(py,b)
                    #ax.scatter(a, b ,s=10,c='red')                                          #Plot Material Points

        px=px[1:(nx*ny)+1]; py=py[1:(nx*ny)+1]
        dVol = np.ones(len(px))*(Vol/(len(px)))                                               #Volume of each material point
        M=np.array((px,py,dVol))                                                            #Array of Material points and dvol
        M=M.transpose()
        return (M)

    if (Pattern == 2):
        w = 1/math.sqrt(3);
        #fig, ax = plt.subplots(1,1,figsize = (5,5))
        px = 0; py=0;

        for i in range (nx):
            #ax.plot([x_0+dx*i,x_0+dx*i],[y_0,l+y_0], c=(0.9, 0.9, 0.9, 1.0) )                             #Plot X lines
            for j in range (ny):        
                #ax.plot([x_0,l+x_0],[y_0+dy*j,y_0+dy*j], c=(0.9, 0.9, 0.9, 1.0) )                         #Plot Y lines
        
                a1 = (x_0+dx*i+x_0+dx*(i+1))/2 - dx*w*0.5; b1 = (y_0+dy*j+y_0+dy*(j+1))/2 - dy*w*0.5;
                if ( math.sqrt(((a1-x_c)**2) + ((b1-y_c)**2)) <= r   ):
                    px = np.append(px,a1)
                    py = np.append(py,b1)
                    #ax.scatter(a1, b1, s=10,c='red')                                                       #Plot Material Points

                a2 = (x_0+dx*i+x_0+dx*(i+1))/2 + dx*w*0.5; b2 = (y_0+dy*j+y_0+dy*(j+1))/2 + dy*w*0.5;
                if ( math.sqrt(((a2-x_c)**2) + ((b2-y_c)**2)) <= r   ):
                    px = np.append(px,a2)
                    py = np.append(py,b2)
                    #ax.scatter(a2, b2, s=10,c='red')                                                       #Plot Material Points
        
                a3 = (x_0+dx*i+x_0+dx*(i+1))/2 - dx*w*0.5; b3 = (y_0+dy*j+y_0+dy*(j+1))/2 + dy*w*0.5;
                if ( math.sqrt(((a3-x_c)**2) + ((b3-y_c)**2)) <= r   ): 
                    px = np.append(px,a3)
                    py = np.append(py,b3)
                    #ax.scatter(a3, b3,s=10,c='red')                                                        #Plot Material Points
        
                a4 = (x_0+dx*i+x_0+dx*(i+1))/2 + dx*w*0.5; b4 = (y_0+dy*j+y_0+dy*(j+1))/2 - dy*w*0.5; 
                if ( math.sqrt(((a4-x_c)**2) + ((b4-y_c)**2)) <= r   ):
                    px = np.append(px,a4)
                    py = np.append(py,b4)
                    #ax.scatter(a4, b4,s=10,c='red')                                                        #Plot Material Points

        px=px[1:(nx*ny*4)+1]; py=py[1:(nx*ny*4)+1]
        dVol = np.ones(len(px))*(Vol/(len(px)))                                                     #Volume of each material point
        M=np.array((px,py,dVol))                                                                 #Array of Material points and dvol
        M=M.transpose()
        return M

    if (Pattern != 1 and Pattern != 2):
        print("Valid entries are only 1 and 2")
        print("*** Nothing returned ***")