#!/usr/bin/python

import matplotlib.pyplot as plt

import csv
import sys
import numpy

def data_from_file(filename):
    file = open(filename, 'rb')
    data=numpy.loadtxt(open(filename,"rb"),delimiter=",",skiprows=0,dtype=int)
    #header=numpy.loadtxt(open(filename,"rb"),delimiter=",",skiprows=2)
    for i, line in enumerate(file):
        if i == 0:
            line0 = line.rstrip('\n')
        elif i == 1:
            line1 = line.rstrip('\n')
        elif i > 1:
            break
    yunit = line0.split(',')[0].rstrip(')').split('(')[1]
    xunit = line0.split(',')[1].split(' ')[2]
    heads = [line.rstrip(' ') for line in line1.split(',')]
    heads[0] = heads[0].lstrip('#')
    #print yunit
    print xunit
    #print heads

    fig=plt.figure()

    n = len(heads)
    x = numpy.array(range(len(data[:,0])))*int(xunit)
        
    for i in range(n):
        ax=fig.add_subplot(n,1,i+1)
        ax.ticklabel_format(useOffset=False)
        ax.set_title(heads[i])
        ax.plot(x,data[:,i])
        ax.set_ylabel(yunit)
        if i+1 is not n:
            ax.get_xaxis().set_ticks([])
        else:
            ax.set_xlabel('s')
    plt.show()
    
        
if __name__=="__main__":
    filename=sys.argv[1]
    data_from_file(filename)
