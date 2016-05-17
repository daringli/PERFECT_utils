import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def diverging_rb_cm():
    cdict_d = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 1.0),
                   (1.0, 0.1, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.1),
                   (0.5, 1.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
    
    cdict = {'red':  ((0.0, 0.0, 0.0),
                      (0.25, 0.0, 0.0),
                      (0.5, 0.8, 1.0),
                      (0.75, 1.0, 1.0),
                      (1.0, 0.4, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (0.25, 0.0, 0.0),
                       (0.5, 0.9, 0.9),
                       (0.75, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),

             'blue':  ((0.0, 0.0, 0.4),
                       (0.25, 1.0, 1.0),
                       (0.5, 1.0, 0.8),
                       (0.75, 0.0, 0.0),
                       (1.0, 0.0, 0.0))
    }
    cdict2 = cdict_d.copy()
    cdict2['alpha'] = ((0.0, 1.0, 1.0),
                #   (0.25,1.0, 1.0),
                   (0.5, 0.3, 0.3),
                #   (0.75,1.0, 1.0),
                   (1.0, 1.0, 1.0))
    plt.register_cmap(name='BlueRedAlpha', data=cdict2)
    blue_red = LinearSegmentedColormap('BlueRedAlpha', cdict2)
    return blue_red
