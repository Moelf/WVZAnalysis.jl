import json
from typing import Dict, List, Optional, Tuple
import numpy as np

import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split


DATA_PATHS = {
    'Signal': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_VVZ.arrow',
    'Other': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_others.arrow',
    'ttZ': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_ttZ.arrow',
    'tWZ': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_tWZ.arrow',
    'tZ': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_tZ.arrow',
    'WZ': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_WZ.arrow',
    'Zgamma': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_Zgamma.arrow',
    'Zjets': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_Zjets.arrow',
    'ZZ': '/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_ZZ.arrow'
}


def load_datasets_from_arrow(
    data_paths: Optional[Dict[str, str]] = None
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Load the datasets from arrow files.

    Parameters
    ----------
    data_paths : Dict[str, str] or `None`
        Dictionary linking source names to their respective `arrow` files.
        Must contain a source labeled `Signal`.
        If `None`, uses default datasets created on 2022-03-01.

    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame]
        Signal and background datasets, in that order.
    '''
    try:
        signal = pd.read_feather(data_paths['Signal'])
        signal['source'] = 'Signal'
        signal['is_signal'] = True
    except KeyError:
        raise KeyError('Data path dictionary must contain a signal source.')

    other_sources = [k for k in data_paths.keys() if k != 'Signal']
    backgrounds = [None] * len(other_sources)
    for i, source in enumerate(other_sources):
        bg = pd.read_feather(data_paths[source])
        bg['source'] = source
        bg['is_signal'] = False

        backgrounds[i] = bg

    background = pd.concat(backgrounds)

    return signal, background


def cut_to_sr(sig: pd.DataFrame, bg: pd.DataFrame,
              sr_index: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Cut the signal and background down to a particular signal region.

    Parameters
    ----------
    sig : pandas.DataFrame
        Signal sample.
    bg : pandas.DataFrame
        Background sample.
    sr_index: int
        Index of the signal region. Absent rescaling, should be
        0 for SF in Z, 1 for SF no Z, and 2 for DF.

    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame]
        Signal and background samples within the desired signal region.
    '''
    bg = bg[bg.SR == sr_index]
    sig = sig[sig.SR == sr_index]

    return sig, bg


def generate_scale_params_file(data_paths: Optional[Dict[str, str]], rescale_features: List[str],
                               json_filepath: str):
    '''
    Given a list of training observables, generate the scale parameters
    for a min/max scaling and write them to a `.json` file.

    Parameters
    ----------
    data_paths : Dict[str, str] or `None`
        Dictionary linking source names to their respective `arrow` files.
        Must contain a source labeled `Signal`.
        If `None`, uses default datasets created on 2022-03-01.
    rescale_features : List[str]
        List of features to rescale.
    json_filepath : str
        Filepath to a `.json` file to which to write the scale parameters.
    '''
    if not json_filepath.endswith('.json'):
        raise ValueError('JSON filepath must be a ".json" file.')

    sig, bg = load_datasets_from_arrow(data_paths)
    full_dataset = pd.concat([sig, bg], ignore_index=True)

    min_max_scaler = preprocessing.MinMaxScaler()
    min_max_scaler.fit(full_dataset[rescale_features])

    scale_params = {}
    scale_params['min'] = {feat: min_max_scaler.min_[i]
                           for i, feat in enumerate(rescale_features)}
    scale_params['scale'] = {feat: min_max_scaler.scale_[i]
                             for i, feat in enumerate(rescale_features)}

    with open(json_filepath, 'w', encoding='utf-8') as json_file:
        json.dump(scale_params, json_file)


def min_max_scale_datasets_from_file(sig: pd.DataFrame,
                                     bg: pd.DataFrame,
                                     feats_to_rescale: List[str],
                                     rescaling_filepath: str):
    '''
    Min/max scale particular features in the signal and background samples
    according to the rescaling parameters listed in a file. Performs rescaling
    in place.

    Parameters
    ----------
    sig : pd.DataFrame
        Signal sample.
    bg : pd.DataFrame
        Background sample.
    feats_to_rescale : List[str]
        List of features to rescale.
    rescaling_filepath : str
        Path of file containing rescaling parameters.
    '''
    with open(rescaling_filepath, 'r', encoding='utf-8') as json_file:
        min_max_scale_params = json.load(json_file)

    for df in [sig, bg]:
        for train_feat in feats_to_rescale:
            scale_val = min_max_scale_params['scale'][train_feat]
            min_val = min_max_scale_params['min'][train_feat]

            df[train_feat] = df[train_feat] * scale_val + min_val


def get_train_test_val_data(data_paths: Dict[str, str],
                            train_feats: List[str],
                            sr_to_train: str,
                            test_prop: float,
                            val_prop: float,
                            rescale_filepath: str,
                            rescale_feats: Optional[List[str]] = None,
                            verbose: bool = True,
                            ) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray,
                                       pd.DataFrame, np.ndarray, np.ndarray,
                                       pd.DataFrame, np.ndarray, np.ndarray]:
    '''
    Load training and validation data, min/max scale training features,
    and split into training, test, and validation sets. Loads data within
    a particular signal region.

    Parameters
    ----------
    data_paths : Dict[str, str]
        Paths of datasets to be loaded, indexed by source names.
        Must contain a source labeled `Signal`.
    train_feats : List[str]
        List of training features to load.
    sr_to_train : str
        Name of the signal region to load. Must be `DF`, `SF_noZ`, `SF_inZ`, or 'full'.
    test_prop : float
        Proportion of events in test sample.
    val_prop : float
        Proportion of events in validation sample.
    rescale_file : str
        Filepath of file containing rescaling parameters.
    rescale_feats : (optional) List[str]
        List of features to rescale. If `None`, all training features will be rescaled.
    verbose : (optional) bool
        Verbosity of output.

    Returns
    -------
    Tuple of DataFrames and ndarrays
        3 groups (train, test, and validation) of 3 samples
        (training features, truth labels, weights).
    '''
    if sr_to_train == 'DF':
        sr_index = 2
    elif sr_to_train == 'SF_noZ':
        sr_index = 1
    elif sr_to_train == 'SF_inZ':
        sr_index = 0
    elif sr_to_train == 'full':
        pass
    else:
        raise KeyError('invalid signal region name')

    sig, bg = load_datasets_from_arrow(data_paths)

    if verbose:
        print('Data loaded...')

    if rescale_feats is not None:
        min_max_scale_datasets_from_file(sig, bg, rescale_feats, rescale_filepath)
    else:
        min_max_scale_datasets_from_file(sig, bg, train_feats, rescale_filepath)

    if verbose:
        print('Data scaled...')

    if sr_to_train != 'full':
        sig, bg = cut_to_sr(sig, bg, sr_index)

        if verbose:
            print('Data cut down to ' + sr_to_train + ' signal region...')

    sig_train, sig_val = train_test_split(sig, test_size=val_prop + test_prop)
    sig_val, sig_test = train_test_split(sig_val, test_size=test_prop / (val_prop + test_prop))

    bg_train, bg_val = train_test_split(bg, test_size=val_prop + test_prop)
    bg_val, bg_test = train_test_split(bg_val, test_size=test_prop / (val_prop + test_prop))

    x_train = pd.concat([sig_train[train_feats], bg_train[train_feats]])
    y_train = np.concatenate([np.ones(len(sig_train)), np.zeros(len(bg_train))])

    n_train_sig = sum(sig_train['wgt'])
    n_train_bg = sum(bg_train['wgt'])

    sig_correction = (n_train_sig + n_train_bg) / (2 * n_train_sig)
    bg_correction = (n_train_sig + n_train_bg) / (2 * n_train_bg)

    w_train = np.concatenate([sig_correction * np.abs(sig_train['wgt']),
                              bg_correction * np.abs(bg_train['wgt'])])

    x_val = pd.concat([sig_val[train_feats], bg_val[train_feats]])
    y_val = np.concatenate([np.ones(len(sig_val)), np.zeros(len(bg_val))])
    w_val = np.concatenate([np.abs(sig_val['wgt']), np.abs(bg_val['wgt'])])

    x_test = pd.concat([sig_test[train_feats], bg_test[train_feats]])
    y_test = np.concatenate([np.ones(len(sig_test)), np.zeros(len(bg_test))])
    w_test = np.concatenate([np.abs(sig_test['wgt']), np.abs(bg_test['wgt'])])

    if verbose:
        print('Splits generated... Finished.')

    return x_train, y_train, w_train, x_test, y_test, w_test, x_val, y_val, w_val
