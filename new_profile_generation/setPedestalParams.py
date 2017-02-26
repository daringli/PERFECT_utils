# specify 3 parameters in the equation below and get the 4th one
# YPed + ddx_YPed * width - YLCFS = 0

tol = 1e-13

def setPedestalParams(YPed=None,ddx_YPed=None, width=None,YLCFS=None):
    noneIndex = [i for i,x in enumerate([YPed,ddx_YPed,width,YLCFS]) if x==None]
    if len(noneIndex) > 1:
        print "setPedestalParams: error: more than one unspecified parameter!"
        raise ValueError('More than one unspecified parameter')
    elif len(noneIndex) == 0:
        # all values are already provided.
        # Check if they are consistent and return them
        if abs(YPed + ddx_YPed * width - YLCFS) <= tol:
            return (YPed,ddx_YPed,width,YLCFS)
        else:
            print "setPedestalParams: error: Inconsistent set of parameters!"
            raise ValueError('Pedestal parameters are not mathematically consistent')
    else:
        noneIndex = noneIndex[0]
        if noneIndex == 0:
            YPed = YLCFS - ddx_YPed * width
        elif noneIndex == 1:
            ddx_YPed = (YLCFS - YPed)/width
        elif noneIndex == 2:
            width = (YLCFS - YPed)/ddx_YPed
        else:
            YLCFS = YPed + ddx_YPed * width
        return (YPed,ddx_YPed,width,YLCFS)

def setPedestalStartStop(XPedStart=None,XPedStop=None,width=None):
    noneIndex = [i for i,x in enumerate([XPedStart,XPedStop,width]) if x==None]
    if len(noneIndex) > 1:
        print "setPedestalParams: error: more than one unspecified parameter!"
        raise ValueError('More than one unspecified parameter')
    elif len(noneIndex) == 0:
        # all values are already provided, so return the inputs
        if abs(XPedStart + width - XPedStop) <= tol:
            return (XPedStart,XPedStop,width)
        else:
            print "setPedestalParams: error: Inconsistent set of parameters!"
            raise ValueError('Pedestal extent are not mathematically consistent')
    else:
        noneIndex = noneIndex[0]
        if noneIndex == 0:
            XPedStart = XPedStop - width
        elif noneIndex == 1:
            XPedStop = XPedStart + width
        else:
            width = XPedStop - XPedStart
        return (XPedStart,XPedStop,width)


if __name__ == "__main__":
    print setPedestalParams(YPed=1,ddx_YPed=-1,width=1) #(1,-1,1,0)
    print setPedestalParams(YPed=1,ddx_YPed=-1,YLCFS=0.5) #(1,-1,0.5,0.5)
    print setPedestalParams(YPed=1.5,width=1,YLCFS=0.5) #(1.5,-1,1,0.5)
    print setPedestalParams(ddx_YPed=-2,width=1,YLCFS=0.5) #(2.5,-2,1,0.5)
    try:
        print setPedestalParams(YPed=1,ddx_YPed=-2,width=1,YLCFS=0.5) # bad
    except ValueError:
        print "inconsistent!"
    print setPedestalStartStop(XPedStart=0,XPedStop=1) #(0,1,1)
    print setPedestalStartStop(XPedStart=0.1,width=0.43) #(0.1,0.53,0.43)
    print setPedestalStartStop(XPedStop=0.7,width=0.23) #(0.47,0.7,0.23)
    try:
        print setPedestalStartStop(XPedStart=0.1,XPedStop=0.7,width=0.23) # bad
    except ValueError:
        print "inconsistent!"
    
