from __future__ import division


import numpy

def diff_matrix(xMin,xMax,N,order=4,periodic=False,endpoint=False):
    ddx=numpy.zeros((N,N))
    if periodic:
        # need to decide whether to include xMax or not
        if endpoint:
            x = [float((xMax-xMin)*i)/(N)+xMin for i in range(1,N+1)]
        else:
            x = [float((xMax-xMin)*i)/(N)+xMin for i in range(0,N)]
        dx=x[1]-x[0] 
    else:
        #includes both end-points
        x = [float((xMax-xMin)*i)/(N-1)+xMin for i in range(0,N)]
    dx=x[1]-x[0]

    if order==4:
        for i in range(2,N-2):
            ddx[i,(i-2)]=1.0/(6*2*dx)
            ddx[i,(i-1)]=-4.0/(3*2*dx)
            ddx[i,(i+1)]=4.0/(3*2*dx)
            ddx[i,(i+2)]=-1.0/(6*2*dx)
            

        if periodic:
            ddx[0, N-2] = 1.0/(6*2*dx)
            ddx[0, N-1] = -4.0/(3*2*dx)
            ddx[0,1]=4.0/(3*2*dx)
            ddx[0,2]=-1.0/(6*2*dx)
            
            ddx[1, N-1] = 1.0/(6*2*dx)
            ddx[1,0]=-4.0/(3*2*dx)
            ddx[1,2]=4.0/(3*2*dx)
            ddx[1,3]=-1.0/(6*2*dx)

            ddx[N-2,N-4]=1.0/(6*2*dx)
            ddx[N-2,N-3]=-4.0/(3*2*dx)
            ddx[N-2,N-1]=4.0/(3*2*dx)
            ddx[N-2, 0] = -1.0/(6*2*dx)

            ddx[N-1,N-3]=1.0/(6*2*dx)
            ddx[N-1,N-2]=-4.0/(3*2*dx)
            ddx[N-1, 0] = 4.0/(3*2*dx)
            ddx[N-1, 1] = -1.0/(6*2*dx)
            
        else:
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
            ddx[i,(i-1)]=-1.0/(2.0*dx)
            ddx[i,(i+1)]=1.0/(2.0*dx)
        if periodic:
            ddx[0,N-1] = -1/(2*dx)
            ddx[0,1] = 1/(2*dx)
            ddx[N-1,0] = 1/(2*dx)
            ddx[N-1,N-2] = -1/(2*dx)
        else:
            ddx[0,0]=-3.0/(2.0*dx)
            ddx[0,1]=2.0/dx
            ddx[0,2]=-1.0/(2.0*dx)

            ddx[-1,-1]=3.0/(2.0*dx)
            ddx[-1,-2]=-2.0/dx
            ddx[-1,-3]=1.0/(2.0*dx)
    else:
        print "ERROR: diff_matrix: Order " + str(order) + " not supported!"
    
    #ddx=numpy.transpose(ddx) 
    return ddx

if __name__=="__main__":
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import matplotlib.pyplot as plt

    Npsi=6
    Ntheta=9
    Nperiods=1
    theta=numpy.linspace(0,Nperiods*2*numpy.pi,Ntheta, endpoint=False)
    psi=numpy.linspace(0,Nperiods*2*numpy.pi,Npsi, endpoint=True)
    print theta[-1]/(numpy.pi)

    ddtheta = diff_matrix(theta[0],theta[0]+Nperiods*2*numpy.pi,N=Ntheta,order=4,periodic=True,endpoint=False)
    
    ddpsi = diff_matrix(psi[0],psi[-1],Npsi,order=4,periodic=False,endpoint=False)
    
    ## 1D #################

    
    z1=numpy.sin(psi)
    z2=numpy.sin(theta)
    
    #dz1dp = numpy.tensordot(ddpsi,z1,([1],[0]))
    #dz2dt = numpy.tensordot(ddtheta,z2,([1],[0]))

    #dz1dp = numpy.dot(ddpsi,z1)
    #dz2dt = numpy.dot(ddtheta,z2)

    dz1dp = numpy.einsum('ij,j',ddpsi,z1)
    dz2dt = numpy.einsum('ij,j',ddtheta,z2)

    #dz1dp=numpy.zeros(Npsi)
    #dz2dt=numpy.zeros(Ntheta)

    #matrix multiplication
    #for i in range(Npsi):
    #    for j in range(Npsi):
    #        dz1dp[i] += ddpsi[i,j]*z1[j]
    #for i in range(Ntheta):
    #    for j in range(Ntheta):
    #        dz2dt[i] += ddtheta[i,j]*z2[j]
    #        print "ddtheta[" + str((i,j)) + "]: " + str( ddtheta[i,j])
    #        print "z[" + str((j)) + "]: " + str(z2[j])
    #    print "dz/dt " + str(dz2dt[i])
    #    print "-------"

    numpy.savetxt("ddtheta.txt",ddtheta)

    analytic_dz1dp=numpy.cos(psi)
    analytic_dz2dt= numpy.cos(theta)


    fig = plt.figure()

    ax1 = fig.add_subplot(2, 1, 1);
    ax2 = fig.add_subplot(2, 1, 2);

    ax1.plot(psi, dz1dp)
    ax1.plot(psi, analytic_dz1dp)

    ax2.plot(theta, dz2dt)
    ax2.plot(theta, analytic_dz2dt)




    plt.show()

    ## 2D #################
    
    P,T=numpy.meshgrid(psi,theta,indexing='ij')
    
    Z = P*numpy.sin(T)
    
    analytic_dZdP = numpy.sin(T)
    analytic_dZdT = P*numpy.cos(T)
    print Z.shape

    dZdP = numpy.einsum('ij,jk',ddpsi,Z)
    dZdT = numpy.einsum('kj,ij',ddtheta,Z)

    print dZdP.shape
    print dZdT.shape

    #plot to verify that analytic and numerical agree well
    fig = plt.figure()

    ax1 = fig.add_subplot(2, 1, 1);
    ax2 = fig.add_subplot(2, 1, 2);
    
    ax1.plot(theta, dZdP[0,:])
    ax1.plot(theta, analytic_dZdP[0,:])
    
    ax2.plot(psi, dZdT[:,Ntheta/2])
    ax2.plot(psi, analytic_dZdT[:,Ntheta/2])

    plt.show()

    ## 3D #################

    Nspecies = 2
    s=numpy.array(range(1,Nspecies+1))
    print s
    P,T,S=numpy.meshgrid(psi,theta,s,indexing='ij')
    print S.shape
    
    Z = P*numpy.sin(T)*S #second species has twice of first
    
    analytic_dZdP = numpy.sin(T)*S
    analytic_dZdT = P*numpy.cos(T)*S

    #dZdP = numpy.tensordot(ddpsi,Z,([1],[0]))
    #dZdT = numpy.transpose(numpy.tensordot(ddtheta,Z,([0],[1])))
    dZdP = numpy.einsum('ij,jkl',ddpsi,Z)
    dZdT = numpy.einsum('kj,ijl',ddtheta,Z)

    print dZdP.shape
    print dZdT.shape

    #plot to verify that analytic and numerical agree well
    fig = plt.figure()

    ax1 = fig.add_subplot(2, 1, 1);
    ax2 = fig.add_subplot(2, 1, 2);
    
    ax1.plot(theta, dZdP[0,:,0])
    ax1.plot(theta, analytic_dZdP[0,:,0])
    
    ax2.plot(psi, dZdT[:,Ntheta/2,1])
    ax2.plot(psi, analytic_dZdT[:,Ntheta/2,1])

    plt.show()
    
