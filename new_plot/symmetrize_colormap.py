from shiftedColorMap import shiftedColorMap
import numpy
import matplotlib
import matplotlib.pyplot as plt
#for __main__
from matplotlib.pyplot import cm
import numpy
from diverging_red_blue_colormap import diverging_rb_cm
from colormap_remove_middle import cm_remove_middle



def invert_cm(cm,name="reversedCM"):
    #inverts a colormap
    new_cm=shiftedColorMap(cm,midpoint=0.5,start=1.0,stop=0.0)
    return new_cm


def symmetrize_colormap(cm,vmin,vmax):
    #if min<0 and max>0, centers the colormap around zero
    # and if abs(min) != abs(max), truncates it to make values with the same absolute value
    # have colors with an equal distance from zero.
    if vmax*vmin<0:
        if vmax>abs(vmin):
            start=(vmax-abs(vmin))/(2.0*vmax) 
            stop=1.0
        else:
            start=0.0
            stop=(abs(vmin)+vmax)/(2.0*abs(vmin)) 
        midpoint=abs(vmin)/float((vmax+abs(vmin)))
        new_cm=shiftedColorMap(cm,midpoint=midpoint,start=start,stop=stop,name="new_map")
    else:
        if vmax>0:
            # truncate lower half of colormap
            start=0.5+0.5*abs(vmin/float(vmax))
            #start = (vmax-abs(vmin))/(2*vmax)
            stop=1.0
        else:
            start=0.0
            stop=0.5-0.5*abs(vmax/float(vmin))
            #stop=(abs(vmin)+vmax)/(2*abs(vmin))
        new_cm=shiftedColorMap(cm,start=start,stop=stop,name="new_map",nomid=True)
    
    return new_cm


if __name__=="__main__":
    #cm=diverging_rb_cm()
    #cm=invert_cm(cm.RdBu)
    cm=invert_cm(cm_remove_middle(cm.RdBu,cut_radius=0.1))
    x=numpy.linspace(0,2*numpy.pi)
    y=numpy.linspace(0,2*numpy.pi)
    X,Y=numpy.meshgrid(x,y)
    z1=1.5*(1+numpy.sin(X-Y))
    z2=2+numpy.sin(X-Y)
    z3=numpy.sin(X-Y)
    cm1=symmetrize_colormap(cm,0,6)
    cm2=symmetrize_colormap(cm,1,3)
    cm3=symmetrize_colormap(cm,-1,1)
    f,(ax1,ax2,ax3)=plt.subplots(3, 1)
    p1=ax1.pcolor(X, Y, z1,rasterized=True,linewidth=0,cmap=cm1)
    p2=ax2.pcolor(X, Y, z2,rasterized=True,linewidth=0,cmap=cm2)
    p3=ax3.pcolor(X, Y, z3,rasterized=True,linewidth=0,cmap=cm3)
    plt.colorbar(p1, ax=ax1)
    plt.colorbar(p2, ax=ax2)
    plt.colorbar(p3, ax=ax3)
    plt.show()
