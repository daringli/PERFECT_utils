from __future__ import division

import matplotlib.pyplot as plt
import numpy as np

def radar_marker_plot(ax,x,y,values,normvalues=None,colors=None,radius=10,markeredgewidth=0.1):
    # based on http://matplotlib.org/examples/api/scatter_piecharts.html
    #generates a radar chart from a value vector and a radius
    #when the value is equal to the norm values, the radius will be equal to radius

    N = len(values)

    if normvalues is not None:
        if len(normvalues) != N:
            raise ValueError("The lenght of normvalues must be equal to the lenght of values!")
    else:
        normvalues = [1]*N

    if colors is not None:
        if type(colors) is not list:
            colors = [colors]
        if len(colors) == 1:
            colors = [colors[0]] * N
        elif len(colors) != N:
            raise ValueError("The lenght of colors must be one equal to the lenght of values!")
    elif N == 4:
        colors =['b','r','g','k']
    else:
        raise ValueError("Automatic colors not supported for N = " + str(N))

    f = 2*np.pi/N # divide unit circle into N parts    

    # to plot

    # shifts can be used to adjust the mid-point of the radar chart
    # relative to the data point.
    # Should be zero unless you know what you are doing
    origin_shift_x = 0
    origin_shift_y = 0
    for i in range(N):        
        marker_x = [origin_shift_x] + np.cos(np.linspace(i*f, (i+1)*f, 10)).tolist()
        marker_y = [origin_shift_y] + np.sin(np.linspace(i*f, (i+1)*f, 10)).tolist()
        marker_xy = list(zip(marker_x, marker_y))
        #s = max(max(marker_x), max(marker_y))
        #size = np.pi*((values[i]/normvalues[i]))**2 * radius
        size = (values[i]/normvalues[i]) * radius

        # plot 
        ax.plot(x, y, marker=(marker_xy, 2),
                   markersize=size, c=colors[i],linestyle='none',markeredgewidth=markeredgewidth)
    

if __name__ == "__main__":

    fig, ax = plt.subplots()
    r=14
    
    x=[0,4,5,0,2,-1]
    y=[0,4,4,2,4,7]
    v=[1,2,3,4]
    nv=[4,4,4,4]    
    marker(ax,x,y,v,nv,radius=r)

    x=[-3,5,6]
    y=[0,-3,-1]
    v=[1,2,4,0]
    nv=[4,4,4,4]
    marker(ax,x,y,v,nv,radius=r)

    x=[-4]
    y=[-4]
    v=[1,2,3,1]
    nv=[4,4,4,4]
    marker(ax,x,y,v,nv,colors='r',radius=r)

    x=[4]
    y=[5]
    v=[2,1,3,1]
    nv=[4,4,4,4]
    marker(ax,x,y,v,nv,colors=['r'],radius=r)

    x=[8]
    y=[8]
    v=[4,3,2,3,4]
    nv=[4,4,4,4,4]
    marker(ax,x,y,v,nv,colors=['r','g','DarkGreen','Gold','Blue'],radius=r)
    
    plt.show()
