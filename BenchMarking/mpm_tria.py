def mpm_tria(centre, dim, Noc, Pattern):

    """
    The function discretises a equilateral/isoceles triangle and returns the material points with available two patterns only.

    Input Parameters:

    Centre: Coordinate of the point were the altutuide meets the base

    dim: Breadth and height of the triangle

    Noc: Number of cells in X and Y direction

    Pattern: Desired choice of dicretisation (Please enter only 1 or 2)
    
    (a material point is created only if the point lies inside the triangle)

    Pattern 1: One material point at the centre of each cell

    Pattern 2: Four material points at each cell distributed based on Guass Legrandge quadrature rule

    Example: MP = mpm_tria([1,2], [3,2], [3,3],2)


    """

    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Circle
    
    #Conversion of data types and assigning variables
    centre = np.array(centre)
    dim = np.array(dim)
    x_c = np.float(centre[0]); y_c = np.float(centre[1]);
    h = np.float(dim[0]); b = np.float(dim[1]);
    nx = np.int(Noc[0]); ny = np.int(Noc[1]);

    x_0 = x_c-(b/2); y_0= y_c;              #Origin of enclosing square
    dx = b/nx; dy = h/ny;                   #Grid spacing in X and Y direction
    Vol = 0.5*b*h;                          #Total Volume

    #function to check if the point is inside triangle
    def tria_check(centre, dim, Coor):
        centre = np.array(centre)
        dim = np.array(dim)
        Coor = np.array(Coor)
    
        p0 = np.zeros(2); p1 = np.zeros(2); p2 = np.zeros(2)
        p0[0] = centre[0] - dim[1]*0.5; p0[1] = centre[1]
        p1[0] = centre[0] + dim[1]*0.5; p1[1] = centre[1]
        p2[0] = centre[0]; p2[1] = centre[1] + dim[0]

        #Barycentric rule
        s = 1/(2*Vol)*(p0[1]*p2[0] - p0[0]*p2[1] + (p2[1] - p0[1])*Coor[0] + (p0[0] - p2[0])*Coor[1]);
        t = 1/(2*Vol)*(p0[0]*p1[1] - p0[1]*p1[0] + (p0[1] - p1[1])*Coor[0] + (p1[0] - p0[0])*Coor[1]);

        if(s>0 and t>0 and 1-s-t>0):
            return 1

    if (Pattern == 1):
        px = 0; py = 0

        for i in range (nx):
            for j in range (ny):        
                a = (x_0+dx*i+x_0+dx*(i+1))/2; b = (y_0+dy*j+y_0+dy*(j+1))/2;
                flag = tria_check(centre, dim, [a,b])
                if ( flag == 1 ):
                    px = np.append(px,a)
                    py = np.append(py,b)
                    
        px=px[1:(nx*ny)+1]; py=py[1:(nx*ny)+1]
        dVol = np.ones(len(px))*(Vol/(len(px)))                                               #Volume of each material point
        M=np.array((px,py,dVol))                                                            #Array of Material points and dvol
        M=M.transpose()
        return (M)

    if (Pattern == 2):
        w = 1/math.sqrt(3);
        px = 0; py=0;

        for i in range (nx):
            for j in range (ny):        
                
                a1 = (x_0+dx*i+x_0+dx*(i+1))/2 - dx*w*0.5; b1 = (y_0+dy*j+y_0+dy*(j+1))/2 - dy*w*0.5;
                flag = tria_check(centre, dim, [a1,b1])
                if ( flag == 1 ):
                    px = np.append(px,a1)
                    py = np.append(py,b1)
                    
                a2 = (x_0+dx*i+x_0+dx*(i+1))/2 + dx*w*0.5; b2 = (y_0+dy*j+y_0+dy*(j+1))/2 + dy*w*0.5;
                flag = tria_check(centre, dim, [a2,b2])
                if ( flag == 1 ):
                    px = np.append(px,a2)
                    py = np.append(py,b2)
                    
                a3 = (x_0+dx*i+x_0+dx*(i+1))/2 - dx*w*0.5; b3 = (y_0+dy*j+y_0+dy*(j+1))/2 + dy*w*0.5;
                flag = tria_check(centre, dim, [a3,b3])
                if ( flag == 1 ): 
                    px = np.append(px,a3)
                    py = np.append(py,b3)
                    
                a4 = (x_0+dx*i+x_0+dx*(i+1))/2 + dx*w*0.5; b4 = (y_0+dy*j+y_0+dy*(j+1))/2 - dy*w*0.5;
                flag = tria_check(centre, dim, [a4,b4])
                if ( flag == 1 ):
                    px = np.append(px,a4)
                    py = np.append(py,b4)
                    
        px=px[1:(nx*ny*4)+1]; py=py[1:(nx*ny*4)+1]
        dVol = np.ones(len(px))*(Vol/(len(px)))                                                     #Volume of each material point
        M=np.array((px,py,dVol))                                                                 #Array of Material points and dvol
        M=M.transpose()
        return M

    if (Pattern != 1 and Pattern != 2):
        print("Valid entries are only 1 and 2")
        print("*** Nothing returned ***")