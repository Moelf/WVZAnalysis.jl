from typing import List
import numpy as np
import pandas as pd

from ..util import util


def left_greedy_scan(data: pd.DataFrame,
                     bin_feature: str,
                     min_count_per_bin: float = 5,
                     n_search_points: int = 100,
                     feat_min: float = 0,
                     feat_max: float = 1,
                     ) -> List[float]:
    '''
    Optimize bin edges by performing a greedy scan from the left.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataset to bin.
    bin_feature : str
        Dataset feature to use for binning (e.g., NN output).
    min_count_per_bin : float (optional, default = 5)
        Minimum number of events (by summed Monte Carlo weight) to keep in each bin.
    n_search_points : int (optional, default = 100)
        Number of cut values to search at each binning step.
    feat_min : float (optional, default = 0)
        Minimum value of the binning feature.
    feat_max : float (optional, default = 1)
        Maximum value of the binning feature.

    Returns
    -------
    List[float]
        List of bin edges.
    '''

    if min_count_per_bin <= 0:
        raise ValueError('Minimum count per bin must be greater than 0.')

    remaining_events = data.copy()
    bin_edges = [feat_min]

    while sum(remaining_events['wgt'] > min_count_per_bin):
        remaining_sigma = 0
        new_edge = feat_max

        for feat_cut in np.linspace(bin_edges[-1], feat_max, n_search_points)[1:]:
            events_to_right = remaining_events[remaining_events[bin_feature] > feat_cut]

            n_sig = np.sum(events_to_right[events_to_right.is_signal == 1]['wgt'])
            n_bg = np.sum(events_to_right[events_to_right.is_signal == 0]['wgt'])

            tmp_sigma = util.significance(n_sig, n_bg)

            current_bin_events = remaining_events[
                (remaining_events[bin_feature] > bin_edges[-1])
                &(remaining_events[remaining_events[bin_feature] < feat_cut])
                ]

            current_bin_count = np.sum(current_bin_events['wgt'])

            if (tmp_sigma > remaining_sigma
                and current_bin_count > min_count_per_bin
                and np.sum(events_to_right['wgt']) > min_count_per_bin):
                remaining_sigma = tmp_sigma
                new_edge = feat_cut

        bin_edges.append(new_edge)
        remaining_events = remaining_events[remaining_events[bin_feature] > new_edge]

    return bin_edges


def right_greedy_scan(data: pd.DataFrame,
                      bin_feature: str,
                      min_count_per_bin: float = 5,
                      n_search_points: int = 100,
                      feat_min: float = 0,
                      feat_max: float = 1,
                      ) -> List[float]:
    '''
    Optimize bin edges by performing a greedy scan from the right.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataset to bin.
    bin_feature : str
        Dataset feature to use for binning (e.g., NN output).
    min_count_per_bin : float (optional, default = 5)
        Minimum number of events (by summed Monte Carlo weight) to keep in each bin.
    n_search_points : int (optional, default = 100)
        Number of cut values to search at each binning step.
    feat_min : float (optional, default = 0)
        Minimum value of the binning feature.
    feat_max : float (optional, default = 1)
        Maximum value of the binning feature.

    Returns
    -------
    List[float]
        List of bin edges.
    '''

    if min_count_per_bin <= 0:
        raise ValueError('Minimum count per bin must be greater than 0.')

    remaining_events = data.copy()
    bin_edges = [feat_max]

    while sum(remaining_events['wgt'] > min_count_per_bin):
        remaining_sigma = 0
        new_edge = feat_min

        for feat_cut in np.linspace(feat_min, bin_edges[0], n_search_points)[:-1]:
            current_bin_events = remaining_events[
                (remaining_events[bin_feature] < bin_edges[0])
                & (remaining_events[bin_feature] > feat_cut)
                ]
            current_bin_count = np.sum(current_bin_events['wgt'])

            n_sig = np.sum(current_bin_events[current_bin_events.is_signal == 1]['wgt'])
            n_bg = np.sum(current_bin_events[current_bin_events.is_signal == 0]['wgt'])

            tmp_sigma = util.significance(n_sig, n_bg)

            if (tmp_sigma > remaining_sigma
                and current_bin_count > min_count_per_bin):
                remaining_sigma = tmp_sigma
                new_edge = feat_cut

        bin_edges.insert(0, new_edge)
        remaining_events = remaining_events[remaining_events[bin_feature] < new_edge]

    return bin_edges
