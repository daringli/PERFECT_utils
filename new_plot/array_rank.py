import numpy

def array_rank(a):
    return len(a.shape)

def arraylist_rank(aol):
    #input may be an array or list of arrays
    if type(aol) is numpy.ndarray:
        return array_rank(aol)
    elif type(aol) is list:
        if len(aol)==0:
            print "arraylist_size: warning: input is empty list"
            return []
        else:
            if all(type(x) is numpy.ndarray for x in aol):
                #all list entries are nd arrays
                return [array_rank(x) for x in aol]
    else:
        print "arraylist_size: error: input is neither list nor array"
        


if __name__=="__main__":
    a=numpy.array([1,2,3])
    print arraylist_rank(a)
    b=numpy.array([[1,2,3],[4,5,6]])
    print arraylist_rank(b)
    c=[a]
    print arraylist_rank(c)
    d=[a,b]
    print arraylist_rank(d)
    e=[]
    print arraylist_rank(e)
    f=1
    print arraylist_rank(f)
