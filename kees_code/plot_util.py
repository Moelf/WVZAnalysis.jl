import os

import matplotlib.pyplot as plt

from atlasify import atlasify


def save_fig(save_dir, save_name):
    if not save_dir.endswith('/'):
        save_dir = save_dir + '/'
        
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    filepath = save_dir + save_name
    plt.savefig(filepath + '.png', pad_inches=0.05, bbox_inches='tight')
    plt.savefig(filepath + '.pdf', pad_inches=0.05, bbox_inches='tight')


def make_nn_output_plot(df_1, df_2, column, df_1_label, df_2_label, weight_col='wgt', density=False, 
                        save=False, save_dir='/', save_name='tmp', title=None, log=True):
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
    
    
def make_nn_output_source_plot(*dfs, column, weight_col='wgt', density=False, 
                               save=False, save_dir='/', save_name='tmp', title=None, log=True, 
                               bins=30):
    vals_to_plot = [df[column] for df in dfs]
    weights_to_plot = [df[weight_col] for df in dfs]
    labels = [df.iloc[0]['source'] for df in dfs]
    
    plt.hist(vals_to_plot, weights=weights_to_plot, bins=bins, density=density, 
             stacked=True, label=labels)

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
    