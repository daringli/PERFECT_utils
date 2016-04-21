import numpy
from array_rank import arraylist_rank

def nstack(array,axis,n):
    rank=arraylist_rank(array)
    if (axis > rank):
        print "nstack: warning: cannot add axis beyond the final position of the input array. Will add an axis to the final position"
        axis=rank
    
    arrays=(numpy.expand_dims(array,axis=axis),)*n
    return numpy.concatenate(arrays,axis=axis)
    

if __name__=="__main__":
    axis=0
    a=numpy.array([1,2,3])
    print nstack(a,axis,2)
    axis=1
    a=numpy.array([1,2,3])
    print nstack(a,axis,2)
    axis=2
    a=numpy.array([1,2,3])
    print nstack(a,axis,2)
    a=numpy.array([[[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]]])
    print a
    b=numpy.array([[1,2,3],[4,5,6]])
    print b
    print nstack(b,1,2)
    print nstack(b,2,2)
