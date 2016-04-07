def update_min_max(self,data,min_max):
    if len(data)==0:
        print "update_min_max: warning, input data is empty."
    newmax=numpy.max(data)
    if newmax>min_max[1]:
        min_max[1]=newmax
    newmin=numpy.min(data)
    if newmin<min_max[0]:
        min_max[0]=newmin
    return min_max
