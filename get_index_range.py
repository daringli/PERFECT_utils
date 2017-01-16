import numpy

def get_index_range(data,drange,ret_range=False,period=None):
    """For increasing data in a numpy.array
    gets the indices that corresponds to the values in ranges.
    The data[returned_indices] will if possible contain the endpoints.
    Period: if specified, the data can be non-increasing if it can be made increasing by adding this multiples of this value each time it decreases"""

    if period is not None:
        data = data.copy() # to not change the original array
        j = 0
        temp = data[0]
        for i in range(1,len(data)):
            if data[i]>temp:
                temp = data[i]
            else:
                j = j+1
            data[i] = data[i] + period*j
        drange =list(drange)
        drange[1] = drange[1] + period*j
        
    if not (numpy.diff(data) >= 0).all():
        print "Error: get_index_range: data is decreasing!"
        if period is not None:
            print "Even though periodicity has been taken into account!"
        return None

    start=numpy.argmin(numpy.fabs(data-drange[0]))
    if drange[0]<data[start]:
        start=max(start-1,0) #cannot go lower than indice zero
    stop=numpy.argmin(numpy.fabs(data-drange[1]))
    if drange[1]>data[stop]:
        stop=min(stop+1,len(data)-1) #cannot go higher than indice len(data)-1
    if ret_range==False:
        return [start,stop]
    else:
        return range(start,stop+1)


if __name__=="__main__":
    data=numpy.linspace(0,1,21)
    print data
    print "this should print 0.9, 0.95, 1:"
    print data[get_index_range(data,[0.9,1],ret_range=True)]
    print data[get_index_range(data,[0.9,255],ret_range=True)]
    print "this should print 0.0, 0.05, 0.1:"
    print data[get_index_range(data,[0,0.1],ret_range=True)]
    print data[get_index_range(data,[-255,0.1],ret_range=True)]

    print "these should print an additional element compared to the above two"
    print data[get_index_range(data,[0.89,1],ret_range=True)]
    print data[get_index_range(data,[0,0.11],ret_range=True)]
    
