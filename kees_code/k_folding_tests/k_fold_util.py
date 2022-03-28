import numpy as np
import matplotlib.pyplot as plt

from atlasify import atlasify


def get_unc_from_bins(data, bin_edges, bin_col):
    bin_unc = [-1] * (len(bin_edges) - 1)
    for i in range(len(bin_edges) - 1):
        data_in_bin = data[(data[bin_col] > bin_edges[i])&(data[bin_col] < bin_edges[i + 1])]
        bin_unc[i] = np.sqrt(np.sum(np.square(data_in_bin['wgt'])))
        
    bin_width = bin_edges[1] - bin_edges[0]
    bin_unc = np.asarray(bin_unc) / (sum(data['wgt']) * bin_width)
    return bin_unc
        

def make_output_overfit_plot(test_sig, test_bg, train_sig, train_bg, nn_out_col):
    # Test BG
    h, b, hist_bg = plt.hist(test_bg[nn_out_col], bins=15, weights=test_bg['wgt'], density=True, 
                             alpha=0.5, label='Test BG')
    bin_centers = (b[1:] + b[:-1]) / 2
    test_bg_unc = get_unc_from_bins(test_bg, b, nn_out_col)
    plt.errorbar(bin_centers, h, test_bg_unc, color=hist_bg[0].get_facecolor(), 
                 alpha=1, linewidth=1, ls='none', capsize=3)
    
    # Test SIG
    h, _, hist_sig = plt.hist(test_sig[nn_out_col], bins=b, weights=test_sig['wgt'], density=True, 
                              alpha=0.5, label='Test SIG')
    test_sig_unc = get_unc_from_bins(test_sig, b, nn_out_col)
    plt.errorbar(bin_centers, h, test_sig_unc, color=hist_sig[0].get_facecolor(), 
                 alpha=1, linewidth=1, ls='none', capsize=3)
    
    # Train BG
    h, _, _ = plt.hist(train_bg[nn_out_col], bins=b, weights=train_bg['wgt'], density=True, 
             histtype='step', color=hist_bg[0].get_facecolor(), alpha=1, label='Train BG')
    train_bg_unc = get_unc_from_bins(train_bg, b, nn_out_col)
    plt.errorbar(bin_centers, h, train_bg_unc, color=hist_bg[0].get_facecolor(), 
                 alpha=1, linewidth=1, ls='none', capsize=3)
    
    # Train SIG
    h, _, _ = plt.hist(train_sig[nn_out_col], bins=b, weights=train_sig['wgt'], density=True, 
             histtype='step', color=hist_sig[0].get_facecolor(), alpha=1, label='Train SIG')
    train_sig_unc = get_unc_from_bins(train_sig, b, nn_out_col)
    plt.errorbar(bin_centers, h, train_sig_unc, color=hist_sig[0].get_facecolor(), 
                 alpha=1, linewidth=1, ls='none', capsize=3)
    
    plt.legend(fontsize=12)
    plt.xlabel('NN output', fontsize=12)
    plt.ylabel('Normalized [a.u.]', fontsize=12)
    atlasify('Internal Simulation', outside=True)