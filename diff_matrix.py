import numpy

def diff_matrix(xMin,xMax,N,order=4):
    ddx=numpy.zeros((N,N))
    x = [float((xMax-xMin)*i)/(N-1)+xMin for i in range(0,N)]
    dx=x[1]-x[0]

    if order==4:
        for i in range(2,N-2):
            ddx[i,i+2]=-1.0/(6*2*dx)
            ddx[i,i+1]=4.0/(3*2*dx)
            ddx[i,i-1]=-4.0/(3*2*dx)
            ddx[i,i-2]=1.0/(6*2*dx)

        ddx[0,0]= -25.0/(12*dx)
        ddx[0,1]= 4.0/(dx)
        ddx[0,2]=-3.0/dx
        ddx[0,3]=4.0/(3*dx)
        ddx[0,4]=-1.0/(4*dx)

        ddx[1,0]= -1.0/(4*dx)
        ddx[1,1]= -5.0/(6*dx)
        ddx[1,2]=3.0/(2*dx)
        ddx[1,3]=-1.0/(2*dx)
        ddx[1,4]=1.0/(12*dx)

        ddx[N-1,N-1]= 25.0/(12*dx)
        ddx[N-1,N-2]= -4.0/(dx)
        ddx[N-1,N-3]=3.0/dx
        ddx[N-1,N-4]=-4.0/(3*dx)
        ddx[N-1,N-5]=1.0/(4*dx)

        ddx[N-2,N-1]= 1.0/(4*dx)
        ddx[N-2,N-2]= 5.0/(6*dx)
        ddx[N-2,N-3]=-3.0/(2*dx)
        ddx[N-2,N-4]=1.0/(2*dx)
        ddx[N-2,N-5]=-1.0/(12*dx)
    elif order==2:
        for i in range(1,N-1):
            ddx[i,i-1]=-1.0/(2.0*dx)
            ddx[i,i+1]=1.0/(2.0*dx)
        ddx[0,0]=-3.0/(2.0*dx)
        ddx[0,1]=2.0/dx
        ddx[0,2]=-1.0/(2.0*dx)
        
        ddx[-1,-1]=3.0/(2.0*dx)
        ddx[-1,-2]=-2.0/dx
        ddx[-1,-3]=1.0/(2.0*dx)
    else:
        print "ERROR: diff_matrix: Order " + str(order) + " not supported!"
        
        
    numpy.transpose(ddx)
    return ddx
