from colour import Color
import numpy

import matplotlib.pyplot as plt

def colortuple(r,g,b):
    tuple255=(r,g,b)
    return Color(rgb=tuple(numpy.array(tuple255)/255.0)).hex

def colorlist(nspecies,nsimul,lg=True):
    if nspecies==3 and nsimul==3:

        
        c11=colortuple(223,101,176)
        c12=colortuple(221,28,119)
        c13=colortuple(189,0,38)

        c21=colortuple(254,153,41)
        c22=colortuple(217,95,14)
        c23=colortuple(153,52,4)

        c31=colortuple(65,182,196)
        c32=colortuple(44,127,184)
        c33=colortuple(37,52,148)     
        

        colorlist=[c11,c12,c13,c21,c22,c23,c31,c32,c33]
    elif nspecies==1 and nsimul==3:
        c1=colortuple(150,150,150)
        c2=colortuple(99,99,99)
        c3=colortuple(37,37,37)
        colorlist=[c1,c2,c3]
    else:
        print "colorlist: warning: unsupported number of species and simulations."
    if lg:
        #colorlist2=[]
        #j=0
        #while j < len(colorlist):
        #    colorlist2=colorlist2+colorlist[j:(nspecies+j)]*2
        #    j=j+nspecies
        #colorlist=colorlist2
        colorlist2=[]
        for c in colorlist:
            colorlist2.append(c)
            colorlist2.append(c)
        colorlist=colorlist2
    #print colorlist
    return colorlist

if __name__=="__main__":
    cl=colorlist(3,3,lg=True)
    for i,c in enumerate(cl):
        plt.plot(i,i,'*',markersize=10,color=c)

    plt.show()
