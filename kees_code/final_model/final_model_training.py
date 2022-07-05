import os
import pickle

import json
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split

from tensorflow.keras.layers import Input, Dense, Dropout
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras import backend as K
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping


def load_sig_bg_datasets():
    sig = pd.read_feather('/home/grabanal/WVZ/gabriel_ML_data/20220301_ELReLMIs54_MUReLMIs31_btag77_VVZ.arrow')
    sig['source'] = 'Signal'
    sig['is_signal'] = True

    bg_others = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                                 + '20220301_ELReLMIs54_MUReLMIs31_btag77_others.arrow'))
    bg_others['source'] = 'Other'

    bg_ttZ = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                              + '20220301_ELReLMIs54_MUReLMIs31_btag77_ttZ.arrow'))
    bg_ttZ['source'] = 'ttZ'

    bg_tWZ = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                              + '20220301_ELReLMIs54_MUReLMIs31_btag77_tWZ.arrow'))
    bg_tWZ['source'] = 'tWZ'

    bg_tZ = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                             + '20220301_ELReLMIs54_MUReLMIs31_btag77_tZ.arrow'))
    bg_tZ['source'] = 'tZ'

    bg_WZ = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                             + '20220301_ELReLMIs54_MUReLMIs31_btag77_WZ.arrow'))
    bg_WZ['source'] = 'WZ'

    bg_Zgamma = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                                 + '20220301_ELReLMIs54_MUReLMIs31_btag77_Zgamma.arrow'))
    bg_Zgamma['source'] = 'Z + gamma'

    bg_Zjets = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                                + '20220301_ELReLMIs54_MUReLMIs31_btag77_Zjets.arrow'))
    bg_Zjets['source'] = 'Z + jets'

    bg_ZZ = pd.read_feather(('/home/grabanal/WVZ/gabriel_ML_data/'
                             + '20220301_ELReLMIs54_MUReLMIs31_btag77_ZZ.arrow'))
    bg_ZZ['source'] = 'ZZ'

    bg = pd.concat([bg_others, bg_ttZ, bg_tWZ, bg_tZ, bg_WZ, bg_Zgamma, bg_Zjets, bg_ZZ])
    bg['is_signal'] = False

    return sig, bg


def rescale_datasets(sig, bg, train_feats):
    with open('/home/kbenkend/WVZ_shared_by_group/WVZAnalysis.jl/kees_code/min_max_scale_params.json') as json_file:
        min_max_scale_params = json.load(json_file)

    for df in [sig, bg]:
        for train_feat in train_feats:
            scale_val = min_max_scale_params['scale'][train_feat]
            min_val = min_max_scale_params['min'][train_feat]

            df[train_feat] = df[train_feat] * scale_val + min_val


def cut_to_sr(sig, bg, sr_index):
    '''
    sr_index: 0 for SF in Z, 0.5 for SF no Z, 1 for DF
    '''
    bg = bg[bg.SR == sr_index]
    sig = sig[sig.SR == sr_index]

    return sig, bg


def get_train_val_data(sr_to_train, val_prop):
    '''
    sr_to_train: 'DF', 'SF_noZ', or 'SF_inZ'
    '''
    if sr_to_train == 'DF':
        sr_index = 1
    elif sr_to_train == 'SF_noZ':
        sr_index = 0.5
    elif sr_to_train == 'SF_inZ':
        sr_index = 0
    else:
        raise KeyError('invalid signal region name')

    sig, bg = load_sig_bg_datasets()
    train_feats = sorted([f for f in sig.columns if f not in ['index', 'wgt', 'is_signal',
                                                              'v_j_btag77', 'v_j_btag60',
                                                              'v_j_btag85', 'v_j_btagCont', 'v_j_btag70',
                                                              'source']])
    print('Using the following training features:')
    print(sorted(train_feats))

    rescale_datasets(sig, bg, train_feats)
    sig, bg = cut_to_sr(sig, bg, sr_index)

    sig_train, sig_val = train_test_split(sig, test_size=val_prop)
    bg_train, bg_val = train_test_split(bg, test_size=val_prop)

    x_train = pd.concat([sig_train[train_feats], bg_train[train_feats]])
    y_train = np.concatenate([np.ones(len(sig_train)), np.zeros(len(bg_train))])

    n_train_sig = sum(sig_train.wgt)
    n_train_bg = sum(bg_train.wgt)

    sig_correction = (n_train_sig + n_train_bg) / (2 * n_train_sig)
    bg_correction = (n_train_sig + n_train_bg) / (2 * n_train_bg)

    w_train = np.concatenate([sig_correction * np.abs(sig_train['wgt']),
                              bg_correction * np.abs(bg_train['wgt'])])

    x_val = pd.concat([sig_val[train_feats], bg_val[train_feats]])
    y_val = np.concatenate([np.ones(len(sig_val)), np.zeros(len(bg_val))])
    w_val = np.concatenate([np.abs(sig_val['wgt']), np.abs(bg_val['wgt'])])

    return x_train, y_train, w_train, x_val, y_val, w_val


def make_model(input_dim, num_nodes, dropout, learn_rate):
    # Generate and fit model
    K.clear_session()
    classifier = Sequential()
    classifier.add(Dense(num_nodes, input_dim=input_dim, activation='relu'))
    classifier.add(Dropout(dropout))
    classifier.add(Dense(num_nodes, activation='relu'))
    classifier.add(Dropout(dropout))
    classifier.add(Dense(num_nodes, activation='relu'))
    classifier.add(Dropout(dropout))
    classifier.add(Dense(1, activation='sigmoid'))

    opt = keras.optimizers.Adam(learning_rate=learn_rate)
    classifier.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])

    return classifier


def train_model(sr_to_train,
                val_prop,
                batch_size=128,
                num_nodes=64,
                dropout=0.1,
                learn_rate=1e-5,
                epochs=10000, patience=200):
    '''
    sr_to_train: 'DF', 'SF_noZ', or 'SF_inZ'
    '''


    x_train, y_train, w_train, x_val, y_val, w_val = get_train_val_data(sr_to_train, val_prop)

    classifier = make_model(input_dim=x_train.shape[1], num_nodes=num_nodes,
                            dropout=dropout, learn_rate=learn_rate)

    es_callback = EarlyStopping(monitor='val_loss', patience=patience,
                                restore_best_weights=True)

    history = classifier.fit(x_train, y_train, sample_weight=w_train,
                             validation_data=(x_val, y_val, w_val),
                             epochs=epochs, batch_size=batch_size,
                             verbose=1, callbacks=[es_callback], shuffle=True)

    # Save model and history
    model_dir = 'models/'
    model_name = 'classifier_' + sr_to_train
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    classifier.save(model_dir + model_name)
    with open(model_dir + model_name + '_history.pkl', 'wb') as file_pi:
        pickle.dump(history.history, file_pi)