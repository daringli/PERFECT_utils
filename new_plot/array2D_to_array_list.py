from array_rank import arraylist_rank
import numpy

def array2D_to_array_list(array,axis=1):
    #array: 2D array to split up into list of arrays
    #axis: axis to split perpendicular to
    other_axis=(axis+1)%2
    if type(arraylist_rank(array)) is list:
        if all(r == 1 for r in arraylist_rank(array)):
            return array
        else:
            print "array2D_to_array_list: warning: array is already list, but arrays have rank>1"
            return array
    elif arraylist_rank(array) == 1:
        return [array]
    elif arraylist_rank(array) == 2:
        return [numpy.take(array,i,axis=axis) for i in range(len(numpy.take(array,0,axis=other_axis)))]
    else:
        print "array2D_to_array_list: warning: arrays of rank>2 not supported"

if __name__ == "__main__":
    a=numpy.array([[1,2,3],[4,5,6]])
    b=array2D_to_array_list(a,0)
    print b
    c=array2D_to_array_list(a,1)
    print c

    d=array2D_to_array_list(b,0)
    print d

    e=array2D_to_array_list(c,0)
    print e

    a2=numpy.array([1,2,3])
    b2=array2D_to_array_list(a2,0)
    print b2
