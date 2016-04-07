import numpy
def boolean_array(array,true_indices):
    zs=numpy.zeros(numpy.shape(array))
    zs[true_indices]=1
    return zs

if __name__ == "__main__":
    a=numpy.array([1,2,3])
    ta=numpy.where(a==1)
    print boolean_array(a,ta)
    b=numpy.array([[1,2,3],[4,5,1]])
    tb=numpy.where(b==1)
    print tb
    print boolean_array(b,tb)

