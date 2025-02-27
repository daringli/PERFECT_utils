from profile import Profile
from generator import Generator
import numpy
from scipy.interpolate import UnivariateSpline
import h5py

import operator

def flat_data(data,axis=0):
    # create an array with axis as first axis,
    # and all other axes flattened as second axis

    Nx_axis = data.shape[axis]
    Nx_other = numpy.prod(data.shape)//Nx_axis

    flat_data = numpy.zeros((Nx_axis,Nx_other))
    for i in range(Nx_axis):
        slices = [slice(None)]*len(data.shape)
        slices[axis]=[i]
        flat_data[i,:] = data[slices].flatten()

    return flat_data

def deflat_data(data,shape,axis=0):
    # invert flat_data
    Nx_axis = data.shape[0]
    shape = list(shape)
    shape[axis] = Nx_axis 
    Nx_other = numpy.prod(shape)/Nx_axis

    deflat_data = numpy.zeros(shape)
    for i in range(Nx_axis):
        for j in range(Nx_other):
            deflat_data[deflat_index(shape,i,j,axis)] = data[i][j]
    return deflat_data

def flat_index(shape,indices,axis=0):
    new_shape = list(shape)
    new_indices = list(indices)
    del new_shape[axis]
    del new_indices[axis]
    # index = (Nn * ...*  N1) *i0 + (Nn * ...*  N2) *i1 +  Nn *i[n-1] + in
    index = sum([numpy.prod(new_shape[(i+1):])*new_indices[i] for i in range(len(new_indices))])
    return index

def deflat_index(shape,lone_index,index,axis=0):
    new_shape = list(shape)
    del new_shape[axis]
    N = len(new_shape)
    indices = [None]*N
    prod = 1
    for i in range(N):
        indices[-1 - i] = (index % (prod * new_shape[-1 - i]))/prod
        index = index - prod * indices[-1 - i]
        prod = prod * new_shape[-1 - i]
    indices.insert(axis,lone_index)
    
    return tuple(indices)

    

class geometryExtrapolatedFromGeometryInput(object):
    """Geometry generator that takes PERFECT input.geometry.h5 files and extrapolates it into the buffer zone"""
    
    def __init__(self,xs,data,quantity=None,interpolation=3,extrapolation=1,xsFirst=(0.0,None),xsLast=(1.0,None),points=None):
        # xs: tuple of position vectors where data is sampled
        # data: sampled data
        # quantity: some quantities have presets
        # interpolation: how to interpolate data
        # extrapolation: how to extrapolate data
        # xLCFS: tuple indicating at which point to extrapolate from
        self.xs = xs
        
        self.quantity = quantity
        if type(data) is numpy.ndarray:
            self.data = data
        elif type(data) is str:
            # interpreting data as a h5file filename
            # quantity species species subset
            self.geometry_filename = data
            self.geometry_file = h5py.File(self.geometry_filename,'r')
            self.group_name =  self.geometry_file["/"].keys()[0]
            self.assert_resolution_from_h5()
            self.data = self.geometry_file["/" + self.group_name + "/" + quantity][()]
            

        if type(interpolation) is int:
            self.interpolation = "smoothingSpline"
            self.order = interpolation
            self.smoothing = 0
        if type(extrapolation) is int:
            self.extrapolation = "polyfit"
            self.order = extrapolation
            self.smoothing = 0

        if points is None:
            self.xsFirst=xsFirst
            self.xsLast=xsLast
        else:
            self.xsLast=tuple([None]*len(xs))
            temp_xsFirst=[None]*len(points)
            for i,p in enumerate(points):
                if p is not None:
                    temp_xsFirst[i] = xs[i][-1-p]
            self.xsFirst = tuple(temp_xsFirst)

    def boundary_index(self,start,step,boundary,op):
        indices=[]
        for i,x in enumerate(self.xs):
            if boundary[i] is not None:
                index = start
                while op(x[index],boundary[i]):
                    index = index + step
                indices.append(index)
            else:
                indices.append(None)
        #indices = indices[-1::-1]
        return indices

    def upper_boundary_index(self):
        step = 1
        start = 0
        boundary = self.xsLast
        op = operator.lt
        boundary_index = self.boundary_index(start,step,boundary,op)
        # for i,bi in enumerate(boundary_index):
        #     if bi is None:
        #         boundary_index[i] = -1
        return boundary_index
            
    def lower_boundary_index(self):
        step = -1
        start = -1
        boundary = self.xsFirst
        op = operator.gt
        boundary_index = self.boundary_index(start,step,boundary,op)
        # for i,bi in enumerate(boundary_index):
        #     if bi is None:
        #         boundary_index[i] = 0
        return boundary_index
        


    def data_inside(self):
        indices2 = self.upper_boundary_index()
        indices1 = [None]*len(indices2)
        return self._data_inside_outside(indices1,indices2)

    def data_outside(self):
        indices1 = self.lower_boundary_index()
        indices2 = [None]*len(indices1)
        return self._data_inside_outside(indices1,indices2)

    def data_inside_outside(self):
        indices1 = self.lower_boundary_index()
        indices2 = self.upper_boundary_index()
        return self._data_inside_outside(indices1,indices2)
        

    def _data_inside_outside(self,indices1,indices2):
        slices = []
        for index1,index2 in zip(indices1,indices2):
            if index1 is not None:
                index1 = index1 + 1
            slices.append(slice(index1,index2,None))
        return self.data[slices]

    def xs_inside_outside(self):
        indices1 = self.lower_boundary_index()
        indices2 = self.upper_boundary_index()
        xs = [None]*len(self.xs)
        for i in range(len(indices1)):
            if indices1[i] is not None:
                xs[i] = self.xs[i][(indices1[i]+1):indices2[i]]
            else:
                xs[i] = self.xs[i][None:indices2[i]]
        return xs
            
    def extrapolate(self,x,axis=0):
        data = self.data_inside_outside()
        
        f_data = flat_data(data,axis)
        xs = self.xs_inside_outside()[axis]
        N = f_data.shape[0]
        N_other = f_data.shape[1]
        results = numpy.zeros((len(x),N_other))
        for i in range(N_other):
            if self.extrapolation == "polyfit":
                p = numpy.polyfit(xs,f_data[:,i],self.order)
                results[:,i] =numpy.polyval(p,x)
        return deflat_data(results,data.shape,axis)
        
    def assert_resolution_from_h5(self):
        Npsi,Ntheta = self.group_name.split('Npsi')[1].split('Ntheta')
        Npsi = int(Npsi)
        Ntheta = int(Ntheta)
        assert (Npsi == len(self.xs[0]))
        assert (Ntheta == len(self.xs[1]))

        
                
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # data
    xp = numpy.linspace(0,1,8) #x points
    yp = numpy.linspace(0,1,6) #y points
    #index order data[x][y]
    data = numpy.transpose(numpy.array([[1,2,3,4,5,6,7,8]]*6)) 
    print xp
    print yp
    print data
    print data[0] # all y values for a given x
    print data[:,0] # all x values for a given y

    xs = (xp,yp)
    xsLast=(0.5,1.0)
    xsFirst=(None,0.2)
    extra=geometryExtrapolatedFromGeometryInput(xs,data,quantity=None,interpolation=3,extrapolation=1,xsFirst=xsFirst,xsLast=xsLast)
    print extra.data_inside()
    print extra.data_outside()
    print extra.data_inside_outside()
    print extra.xs_inside_outside()


    print "##########"

    shape = (3,2,4)
    axis = 0
    a=numpy.ones(shape)*numpy.arange(1,1+shape[2])

    # check that it works to flatten data
    # and that we know how to map flat and nonflat indices
    f_data = flat_data(a,axis)
    df_data = deflat_data(f_data,shape,axis)
    for i0 in range(shape[0]):
        for i1 in range(shape[1]):
            for i2 in range(shape[2]):
                indices = i0,i1,i2
                index = flat_index(shape,indices,axis)
                assert(f_data[indices[axis],index] == a[indices[0],indices[1],indices[2]])
                deflat = deflat_index(shape,indices[axis],index,axis)
                assert(deflat == indices )
                assert((a == df_data).all())

    print "##########"

    x_extra = numpy.linspace(0,2,8)
    print extra.extrapolate(x_extra,0)

    print "##########"
    extra=geometryExtrapolatedFromGeometryInput(xs,data,quantity=None,interpolation=3,extrapolation=1,points = (3,None))
    print extra.data_inside_outside()
    print extra.xs_inside_outside()

