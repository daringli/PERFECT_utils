from colour import Color
import numpy

import matplotlib.pyplot as plt

def colortuple(r,g,b):
    tuple255=(r,g,b)
    return Color(rgb=tuple(numpy.array(tuple255)/255.0)).hex

def colorlist(clist,lg=True):
    #recommended colors:
    #c11=colortuple(223,101,176)
    #c12=colortuple(221,28,119)
    #c13=colortuple(189,0,38)

    #c21=colortuple(254,153,41)
    #c22=colortuple(217,95,14)
    #c23=colortuple(153,52,4)

    #c31=colortuple(65,182,196)
    #c32=colortuple(44,127,184)
    #c33=colortuple(37,52,148) 
    
    if lg:
        colorlist=[]
        for c in clist:
            colorlist.append(c)
            colorlist.append(c)
        clist=colorlist
    return clist

if __name__=="__main__":
    cl=colorlist(3,3,lg=True)
    for i,c in enumerate(cl):
        plt.plot(i,i,'*',markersize=10,color=c)

    plt.show()
