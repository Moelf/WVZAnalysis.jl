import os
import pickle
from typing import Tuple
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend as K
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.models import Sequential
from tensorflow.keras.callbacks import EarlyStopping

physical_devices = tf.config.list_physical_devices('GPU') 
tf.config.experimental.set_memory_growth(physical_devices[0], True)


def make_model(input_dim: int, num_nodes: int, dropout: float, learn_rate: float) -> Sequential:
    '''
    Generate and compile a 3-layer Keras model using the Adam optimizer and binary
    cross-entropy loss.

    Parameters
    ----------
    input_dim : int
        Number of input features
    num_noes : int
        Number of nodes per layer
    dropout : float
        Dropout per layer
    learn_rate : float
        Learning rate for the Adam optimizer

    Returns
    -------
    keras.Sequential
        Compiled model
    '''
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


def make_and_train_model(training_data: Tuple[pd.DataFrame, np.ndarray, np.ndarray],
                         validation_data: Tuple[pd.DataFrame, np.ndarray, np.ndarray],
                         batch_size: int = 128,
                         num_nodes: int = 64,
                         dropout: float = 0.1,
                         learn_rate: float = 1e-5,
                         epochs: int = 10000,
                         patience: int = 200,
                         model_dir: str = 'models/',
                         model_name: str = 'classifier'):
    '''
    Generate a 3-layer model and train it on a particular dataset. Saves the model to disk and the
    training history as a `.pkl` file.

    Parameters
    ----------
    training_data : (pandas.DataFrame, numpy.ndarray, numpy.ndarray)
        Training data in the form (training features, training labels, training weights).
    validation_data : (pandas.DataFrame, numpy.ndarray, numpy.ndarray)
        Validation data in the form (training features, training labels, training weights).
    batch_size: int (default 128)
        Batch size for training.
    num_nodes : int (default 64)
        Number of nodes per layer.
    dropout : float (default 0.1)
        Dropout per layer.
    learn_rate : float (default 1e-5)
        Learning rate for the Adam optimizer.
    epochs : int (default 10000)
        Number of training epochs.
    patience : int (default 200)
        Patience on early stopping during training.
    model_dir : str (default `models/`)
        Directory to save the model.
    model_name : str (default `classifier`)
        Name to which to save the model weights and training history.
    '''


    x_train, y_train, w_train = training_data
    x_val, y_val, w_val = validation_data

    classifier = make_model(input_dim=x_train.shape[1], num_nodes=num_nodes,
                            dropout=dropout, learn_rate=learn_rate)

    es_callback = EarlyStopping(monitor='val_loss', patience=patience,
                                restore_best_weights=True)

    history = classifier.fit(x_train, y_train, sample_weight=w_train,
                             validation_data=(x_val, y_val, w_val),
                             epochs=epochs, batch_size=batch_size,
                             verbose=1, callbacks=[es_callback], shuffle=True)

    # Save model and history
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    classifier.save(model_dir + model_name)
    with open(model_dir + model_name + '_history.pkl', 'wb') as file_pi:
        pickle.dump(history.history, file_pi)