#a naive peakfinder without filters
import numpy
from get_index_range import get_index_range
from diff_matrix import diff_matrix
import matplotlib.pyplot as plt

def PPV(data,x,rangex,tol=1e-10,order=4,mode="default",retmode="y"):
    (peak_y,peak_x,peak_i) = peakfinder(data, x, rangex, tol, order,mode="peaks")
    (dip_y,dip_x,dip_i) = peakfinder(data, x, rangex, tol, order,mode="dips")
    y = numpy.concatenate((peak_y,dip_y))
    oldx=x
    x = numpy.concatenate((peak_x,dip_x))
    indices = numpy.concatenate((peak_i,dip_i))
    # if no peaks or dips are found
    
    this_y=None
    if len(y) == 0:
        print "WARNING: NO PEAK OR DIP WAS FOUND IN xrange: " + str(rangex)
        return float("nan")

    if mode == "default":
        # find most positive peak or most negative dip
        this_y=numpy.concatenate((peak_y,-dip_y))
    elif mode == "reverse":
        # find most negative peak or most positive dip
        # (this one is weird)
        this_y=numpy.concatenate((-peak_y,dip_y))
    elif mode == "positive":
        # find most positive peak
        if len(peak_y) >0:
            this_y=numpy.concatenate((peak_y,0*dip_y))
        else:
            return float("nan")
    elif mode == "negative":
        # find most negative dip
        if len(dip_y) >0:
            this_y=numpy.concatenate((0*peak_y,-dip_y))
        else:
            return float("nan")
    else:
        raise ValueError("Invalid mode: " + retmode)
        
    i=numpy.argmax(this_y)
    if retmode == "y":
        return y[i]
    if retmode == "absy":
        return numpy.fabs(y[i])
    elif retmode == "x":
        return x[i]
    elif retmode == "i":
        return indices[i]
    else:
        raise ValueError("Invalid retmode: " + retmode)

def PPV_x(data,x,rangex,tol=1e-10,order=4,mode="default",retmode="x"):
    # wrapper around PPV with retmode="x". For backwards compatiblity.
    return PPV(data,x,rangex,tol,order,mode,retmode)


def peakfinder(data, x, rangex, tol=1e-10,order=4,mode="both"):
    Nx = len(x)
    ddx=diff_matrix(x[0],x[-1],Nx,order)
    ddx_data = numpy.tensordot(ddx,data,([1],[0]))
    [start,stop] = get_index_range(x,rangex)
    (maxes,mines) = change_sign(ddx_data[start:stop],tol)
    if mode == "both":
        ret = (numpy.concatenate((data[start:stop][maxes],data[start:stop][mines])),numpy.concatenate((x[start:stop][maxes],x[start:stop][mines])),numpy.concatenate((maxes,mines)))        
    elif mode == "peaks":
        ret = (data[start:stop][maxes],x[start:stop][maxes],maxes)
    elif mode == "dips":
        ret = (data[start:stop][mines],x[start:stop][mines],mines)
    else:
        raise ValueError("Invalid mode: " + mode)
    return ret
    
def change_sign(data,tol=1e-10):
    maxes = []
    mines = []
    past_state = None
    for i,d in enumerate(data):
        if d > tol:
            state = "rising"
        elif d < -tol:
            state = "falling"
        else:
            state = "constant"
        #print str(i) + ": " + state
        if past_state is not None:
            if past_state == "rising" and state == "falling":
                maxes.append(i)
            elif past_state == "falling" and state == "rising":
                mines.append(i)
        if state != "constant":
            past_state = state
    return (maxes,mines)

def sweet_visualz(x,y):
    # testz peakfinder and produces sweet visuals using matplotlib(tm)
    rangex = [x[0],x[-1]]
    (y_peak,x_peak,i_peak) = peakfinder(y, x, rangex)
    print "peak/dip locations: " + repr(x_peak)
    print "peak/dip values: " + repr(y_peak)
    ax=plt.gca()
    p1=ax.plot(x,y)
    ax.plot(x_peak,y_peak,'o',color=p1[-1].get_color())

if __name__ == "__main__":
    data = numpy.array([0, 1, 2, 3, 2, 4, 1, 0, -1, -2, 0, 2,-10,-11,-10,-9])
    x = numpy.array([-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,14])
    (mines,maxes) = change_sign(data)
    print mines # [11]
    print maxes # [8]

    # new
    rangex = [0,14]
    maxes = peakfinder(data, x, rangex, 1e-10,4)
    
    print maxes # [3,4]
    print PPV(data, x, rangex, 1e-10,4,mode="default") # 4

    sweet_visualz(x,data)
    plt.show()
