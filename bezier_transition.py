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

    
def bezier_transition(funcListIn,pointList,pairList,x):
    #funcList: list of 1D functions to transition between
    #pointList: list of points that define the transition between regions, and control points of Bezier curves
    #pairList: list of pair of x distances that describe where the Bezier transition starts

    #calculate breakpoints from the point and pair list
    funcList=list(funcListIn)
    breakpoints=[-float("inf")]
    for point,pair in zip(pointList,pairList):
        breakpoints.append(point-pair[0])
        breakpoints.append(point+pair[1])
    breakpoints.append(float("inf"))

    #add transition regions to the funclist
    k=0 #since we will add elements to funclist
    for i in range(len(pointList)):
        x0i=pointList[i]-pairList[i][0]
        x1i=pointList[i]
        x2i=pointList[i]+pairList[i][1]
        
        P0=numpy.array([x0i,funcList[k+i](x0i)])
        P1=numpy.array([x1i,funcList[k+i](x1i)])
        P2=numpy.array([x2i,funcList[k+i+1](x2i)])
   
        def B(P0,P1,P2,x0i,x2i):
            return lambda x: bezier_func_y(P0,P1,P2)((x-x0i)/(x2i-x0i))
        
        funcList.insert(1+2*i, B(P0,P1,P2,x0i,x2i))
        k=k+1

    #create one big function from the final funclist
   
    bigf = []
    bigf.append(lambda x: 0)
    for i in range(len(funcList)):
        bigf.append(lambda x,i=i : bigf[i](x) + characteristic(breakpoints[i],breakpoints[i+1],x)*funcList[i](x))
    return bigf[-1](x)

def derivative_bezier_transition(funcListIn,derivListIn,pointList,pairList,x):
    # RETURNS BOTH FUNCTION AND DERIVATIVE
    # SINCE FUNCTION WAS NEEDED TO CALCULATE THE TRANSITION CURVES
    
    #funcList: list of 1D functions to transition between
    #derivList: list of 1D their derivatives
    #pointList: list of points that define the transition between regions, and control points of Bezier curves
    #deltaPairs: list of pair of x distances that describe where the Bezier transition starts
    #print funcList
    funcList=list(funcListIn)
    derivList=list(derivListIn)
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
        def B(P0,P1,P2,x0i,x2i):
            return lambda x: ddt_bezier_func_y(P0,P1,P2)((x-x0i)/(x2i-x0i))*(1/(x2i-x0i))
        def C(P0,P1,P2,x0i,x2i):
            return lambda x: bezier_func_y(P0,P1,P2)((x-x0i)/(x2i-x0i))
        
        derivList.insert(1+2*i, B(P0,P1,P2,x0i,x2i))
        #need the function to actually calculate the next y point
        funcList.insert(1+2*i, C(P0,P1,P2,x0i,x2i))
        k=k+1

    bigf = []
    bigdfdx = []
    bigf.append(lambda x: 0)
    bigdfdx.append(lambda x: 0)
    for i in range(len(funcList)):
        bigf.append(lambda x,i=i : bigf[i](x) + characteristic(breakpoints[i],breakpoints[i+1],x)*funcList[i](x))
        bigdfdx.append(lambda x,i=i : bigdfdx[i](x) + characteristic(breakpoints[i],breakpoints[i+1],x)*derivList[i](x))

    #2016-08-03: now returns the function, rather than the function sampled at x
    return (bigf[-1],bigdfdx[-1])

        
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
    
