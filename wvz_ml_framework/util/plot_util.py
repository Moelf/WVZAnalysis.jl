import os
from typing import List, Optional, Union

import matplotlib.pyplot as plt

from atlasify import atlasify


COLOR_SCHEME = {'WVZ': '#ff0000', #red
                'ZZ': '#32ffff', #light blue
                'Zjets': '#03ff00', #light green
                'ttZ': '#5ad454', #dark green
                'Zgamma': '#2900ff', #dark blue
                'others': '#ff00ff', #pink
                'tZ': '#660066', #purple
                'WZ': '#cf5f61', #orange
                'tWZ': '#cccccc', #grey
               }


def save_fig(save_dir: str, save_name):
    '''
    Save a figure in a particular location. Automatically applies the correct padding.
    Saves a `.pdf` and `.png` file for each image.

    Save directory and name are separated to ensure that the directory is
    generated correctly.

    Parameters
    ----------
    save_dir : str
        Directory in which to save the figure.
    save_dir : str
        Filename for the figure.
    '''
    if not save_dir.endswith('/'):
        save_dir = save_dir + '/'

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    filepath = save_dir + save_name
    plt.savefig(filepath + '.png', pad_inches=0.05, bbox_inches='tight')
    plt.savefig(filepath + '.pdf', pad_inches=0.05, bbox_inches='tight')


def make_nn_output_plot(df_1, df_2, column, df_1_label, df_2_label, weight_col='wgt', density=False,
                        save=False, save_dir='/', save_name='tmp', title=None, log=True):
    '''
    Plot the output of a neural network.
    '''
    _, b, _ = plt.hist(df_1[column], weights=df_1[weight_col],
                       bins=30, alpha=0.5, label=df_1_label, density=density)
    plt.hist(df_2[column], weights=df_2[weight_col], bins=b, alpha=0.5,
             label=df_2_label, density=density)

    if log:
        plt.yscale('log')

    plt.xlabel('NN Output', fontsize=12)

    if density:
        plt.ylabel('Normalized [a.u.]', fontsize=12)
    else:
        plt.ylabel('Events', fontsize=12)

    if title is not None:
        plt.title(title, fontsize=14, loc='right')

    plt.legend(fontsize=12)

    atlasify('Internal Simulation', outside=True)
    plt.minorticks_on()

    if save:
        save_fig(save_dir, save_name)

    plt.show()


def make_nn_output_source_plot(*dfs, colors: Optional[List[str]], column: str,
                               weight_col: str ='wgt', density: bool = False,
                               save: bool = False,
                               save_dir: str = '/', save_name: str = 'tmp',
                               title: Optional[str] = None,
                               log: bool = True,
                               bins: Union[int, List[float]] = 30):
    '''
    Plot the output of a neural network with the different background components separated
    and colored. This is the "typical" stacked histogram plot that we are used to.

    Parameters
    ----------
    dfs
        Dataframes containing each background source. Different sources must be split into
        different dataframes. The dataframe must have a `source` column listing the background
        source; this will be used for the plot label.
    colors : List[str] or `None`
        List of color values in order of the background source dataframes.
        If `None` is passed, uses a default colorscheme based on the Run 1
        triboson analysis, but the order might be wrong.
    column : str
        Column (i.e., variable) of the dataframe to plot.
    weight_col : str
        Column of the dataframe containing event weights.
    density : (optional) bool, default `False`
        If `True`, plots the density of events (i.e., the normalized event count
        distribution). If `False`, plots event counts.
    save : (optional) bool, default `False`
        Whether to save plot to file.
    save_dir : (optional) str, default `'/'`
        Directory to which to save output plot.
    save_name : (optional) str, default `'tmp'`
        Filename for saved output plot.
    title : (optional) str or `None`, default `None`
        Title for the plot.
    log : (optional) bool, default `True`
        Whether to log-scale the vertical axis.
    bins : (optional) int or List[float], default `30`
        Number of bins to plot or list of bin edges.
    '''
    if colors is None:
        colors = [COLOR_SCHEME[k] for k in COLOR_SCHEME.keys()]

    vals_to_plot = [df[column] for df in dfs]
    weights_to_plot = [df[weight_col] for df in dfs]
    labels = [df.iloc[0]['source'] for df in dfs]

    plt.hist(vals_to_plot, weights=weights_to_plot, bins=bins, density=density,
             stacked=True, label=labels, color=colors, histtype='stepfilled', edgecolor='k')

    if log:
        plt.yscale('log')

    plt.xlabel('NN Output', fontsize=12)

    if density:
        plt.ylabel('Normalized [a.u.]', fontsize=12)
    else:
        plt.ylabel('Events', fontsize=12)

    if title is not None:
        plt.title(title, fontsize=14, loc='right')

    plt.legend(fontsize=12, bbox_to_anchor=(1, 1))

    atlasify('Internal Simulation', outside=True)
    plt.minorticks_on()

    if save:
        save_fig(save_dir, save_name)

    plt.show()
