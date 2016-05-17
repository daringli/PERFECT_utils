import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib
import numpy as np

def cm_remove_middle(cm,midpoint=0.5,cut_radius=0.2,name="my_cm"):
    #removes the middle of a colormap
    # so that region centered around 0 may have sharper contrast

    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
    donext=False
    reg_index = np.linspace(0, 1, 100)
    dri=reg_index[1]-reg_index[0]
    slope1=1/(1-2*cut_radius)
    slope2=0.5/(0.5-cut_radius)
    offset=1-slope2
    for ri in reg_index:
        if (ri > midpoint+cut_radius) or (ri < midpoint - cut_radius):
            r, g, b, a = cm(ri)
            if ri < midpoint :
                cdict['red'].append((ri*slope1, r, r))
                cdict['green'].append((ri*slope1, g, g))
                cdict['blue'].append((ri*slope1, b, b))
                cdict['alpha'].append((ri*slope1, a, a))
                #print str(ri) +":" + str(ri*slope1)
            else:
                cdict['red'].append((offset+ri*slope2, r, r))
                cdict['green'].append((offset+ri*slope2, g, g))
                cdict['blue'].append((offset+ri*slope2, b, b))
                cdict['alpha'].append((offset+ri*slope2, a, a))
                #print str(ri) +":" + str(offset+ri*slope2)
        if ri < (midpoint+dri) and ri > (midpoint-dri):
            r2, g2, b2, a2 = cm(midpoint)
            cdict['red'].append((midpoint, r2, r2))
            cdict['green'].append((midpoint, g2, g2))
            cdict['blue'].append((midpoint, b2, b2))
            cdict['alpha'].append((midpoint, a2, a2))
        
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict,N=512)
    plt.register_cmap(cmap=newcmap)
    return newcmap
    
if __name__=="__main__":
    cm=cm.RdBu
    
    x=np.linspace(0,1)
    y=np.linspace(0,1)
    X,Y=np.meshgrid(x,y)

    Z=X+Y
    plt.subplot(2,1,1)
    plt.pcolor(X,Y,Z,cmap=cm)
    plt.title("RdBu")
    plt.subplot(2,1,2)
    plt.pcolor(X,Y,Z,cmap=remove_middle(cm))
    plt.title("RdBu middle cut")


    plt.show()
