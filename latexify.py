import numpy
import matplotlib

def latexify(fig_width=None, fig_height=None, columns=1):
    # code from http://nipunbatra.github.io/2014/08/latexify/
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    
    if fig_width is None:
        assert(columns in [1,2])
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (numpy.sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    # MAX_HEIGHT_INCHES = 8.0
    # if fig_height > MAX_HEIGHT_INCHES:
    #     print("WARNING: fig_height too large:" + str(fig_height) + 
    #           "so will reduce to" + str(MAX_HEIGHT_INCHES) + "inches.")
    #     fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'pdf',
              #'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 8, # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              #'text.fontsize': 8, # was 10
              'text.color': 'black',
              'legend.fontsize': 8, # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              #'figure.max_num_figures' : 40,
              #'font.family': 'serif',
              #'font.serif': ["Latin modern"]
              'savefig.dpi': 200
    }

    matplotlib.rcParams.update(params)
