import numpy


def bezier_func_y(P0,P1,P2):
    def By(t): return (1-t)*((1-t)*P0[1]+t*P1[1])+t*((1-t)*P1[1]+t*P2[1])
    return By

def ddt_bezier_func_y(P0,P1,P2):
    def By(t): return 2*(t-1)*P0[1] + 2*(1-2*t)*P1[1] + 2*t*P2[1]
    return By


def bezier_func_x(P0,P1,P2):
    def Bx(t): return (1-t)*((1-t)*P0[0]+t*P1[0])+t*((1-t)*P1[0]+t*P2[0])
    return Bx


def step(x):
    #from http://stackoverflow.com/questions/15121048/does-a-heaviside-step-function-exist, works with boty numpy arrays and numbers
    return 1 * (x > 0)

def characteristic(a,b,x):
    return step(x-a)-step(x-b)

    
def bezier_transition(funcList,pointList,pairList,x):
    #funcList: list of 1D functions to transition between
    #pointList: list of points that define the transition between regions, and control points of Bezier curves
    #deltaPairs: list of pair of x distances that describe where the Bezier transition starts
    #print funcList
    breakpoints=[-float("inf")]
    for point,pair in zip(pointList,pairList):
        breakpoints.append(point-pair[0])
        breakpoints.append(point+pair[1])
    breakpoints.append(float("inf"))
    
    k=0 #since we will add elements to funclist
    for i in range(len(pointList)):
        x0i=pointList[i]-pairList[i][0]
        #print x0i
        x1i=pointList[i]
        #print x1i
        x2i=pointList[i]+pairList[i][1]
        #print x2i
        #print "In Bezier trans"
        #print i+k
        #print funcList[k+i]
        #print x0i
        #print funcList[k+i](0.99)
        
        P0=numpy.array([x0i,funcList[k+i](x0i)])
        P1=numpy.array([x1i,funcList[k+i](x1i)])
        P2=numpy.array([x2i,funcList[k+i+1](x2i)])
        #plt.plot(P0[0],P0[1],'o')
        #plt.plot(P1[0],P1[1],'o')
        #plt.plot(P2[0],P2[1],'o')
        def B(P0,P1,P2,x0i,x2i):
            return lambda x: bezier_func_y(P0,P1,P2)((x-x0i)/(x2i-x0i))
        #xtemp=numpy.linspace(x0i,x2i)
        #plt.plot(x,B(P0,P1,P2,x0i,x2i)(x))
        funcList.insert(1+2*i, B(P0,P1,P2,x0i,x2i))
        k=k+1

    fvalue=[0]*len(x)
    for i in range(len(funcList)):
        def f(p1,p2,i,x):
            return characteristic(p1,p2,x)*funcList[i](x)
        #print "Breakpoints:"
        #print breakpoints[i]
        #print breakpoints[i+1]
        #print funcList[i]
        #print ""
        fvalue=fvalue+f(breakpoints[i],breakpoints[i+1],i,x)

    return fvalue

def derivative_bezier_transition(funcList,derivList,pointList,pairList,x):
    #funcList: list of 1D functions to transition between
    #derivList: list of 1D their derivatives
    #pointList: list of points that define the transition between regions, and control points of Bezier curves
    #deltaPairs: list of pair of x distances that describe where the Bezier transition starts
    #print funcList
    breakpoints=[-float("inf")]
    for point,pair in zip(pointList,pairList):
        breakpoints.append(point-pair[0])
        breakpoints.append(point+pair[1])
    breakpoints.append(float("inf"))
    
    k=0 #since we will add elements to funclist
    for i in range(len(pointList)):
        x0i=pointList[i]-pairList[i][0]
        x1i=pointList[i]
        x2i=pointList[i]+pairList[i][1]
        
        P0=numpy.array([x0i,funcList[k+i](x0i)])
        P1=numpy.array([x1i,funcList[k+i](x1i)])
        P2=numpy.array([x2i,funcList[k+i+1](x2i)])
        #plt.plot(P0[0],P0[1],'o')
        #plt.plot(P1[0],P1[1],'o')
        #plt.plot(P2[0],P2[1],'o')
        def B(P0,P1,P2,x0i,x2i):
            return lambda x: ddt_bezier_func_y(P0,P1,P2)((x-x0i)/(x2i-x0i))*(1/(x2i-x0i))
        #xtemp=numpy.linspace(x0i,x2i)
        #plt.plot(x,B(P0,P1,P2,x0i,x2i)(x))
        derivList.insert(1+2*i, B(P0,P1,P2,x0i,x2i))
        k=k+1

    fvalue=[0]*len(x)
    for i in range(len(funcList)):
        def f(p1,p2,i,x):
            return characteristic(p1,p2,x)*derivList[i](x)
        fvalue=fvalue+f(breakpoints[i],breakpoints[i+1],i,x)

    return fvalue


        
if __name__=='__main__':
    #for sanity checking
    import matplotlib.pyplot as plt

    flist=[lambda x: 15-x, lambda x:10-2*(x-5),lambda x:0-(x-10)]
    dlist=[[1,1],[1,1]]
    plist=[5,10]
    #flist=[lambda x: 15-x, lambda x:10-2*(x-5)]
    #dlist=[[1,1]]
    #plist=[5]

    plotNum = 1

    numRows = 2
    numCols = 1
    fig = plt.figure()
    ax = fig.add_subplot(numRows, numCols, plotNum); plotNum += 1
    
    x=numpy.linspace(-5,20,100)
    y=bezier_transition(flist,plist,dlist,x)
    plt.plot(x,y)
    """
    ax = fig.add_subplot(numRows, numCols, plotNum); plotNum += 1
    y=step(x)
    plt.plot(x,y)

    ax = fig.add_subplot(numRows, numCols, plotNum); plotNum += 1
    y=characteristic(3,5,x)
    plt.plot(x,y)
"""
    ax = fig.add_subplot(numRows, numCols, plotNum); plotNum += 1

    t=numpy.linspace(0,1)
    P0=numpy.array([9,0])
    P1=numpy.array([10,-2])
    P2=numpy.array([13,-3])
    y=bezier_func_y(P0,P1,P2)(t)
    x=bezier_func_x(P0,P1,P2)(t)
    plt.plot(x,y)    
    """
    ax = fig.add_subplot(numRows, numCols, plotNum); plotNum += 1
    
    plt.plot(P0,'o') #does not plot a single point!
    plt.plot(P1[0],P1[1],'o')    
    plt.plot(P2[0],P2[1],'o')    
"""
    plt.show()
    
