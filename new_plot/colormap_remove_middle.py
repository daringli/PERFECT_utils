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
    reg_index = np.linspace(0, 1, 257)
    for ri in reg_index:
        if (ri > midpoint+cut_radius) or (ri < midpoint - cut_radius):
            #print ri
            r, g, b, a = cm(ri)
            cdict['red'].append((ri, r, r))
            cdict['green'].append((ri, g, g))
            cdict['blue'].append((ri, b, b))
            cdict['alpha'].append((ri, a, a))
        elif ri==midpoint:
            r, g, b, a = cm(ri)
            cdict['red'].append((ri, r, r))
            cdict['green'].append((ri, g, g))
            cdict['blue'].append((ri, b, b))
            cdict['alpha'].append((ri, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
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
