"""
Tamar Ervin
Date: July 22, 2021

Neural Network outline for looking at shape changes
in residual CCFs
- the residual CCFs should be normalized
    - (ccf - np.median(ccf)) / np.std(ccf)
- also current issue of the residuals looking funky i think
"""

import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

import tensorflow as tf
import tensorflow.keras.backend as BK
from tensorflow.keras.layers import *
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau

from tamar.tools.settings import CsvDir, Config
import matplotlib.pyplot as plt


# tensorflow dataset functions
def load_image_train(datapoint):
    input_image = tf.reshape(datapoint, (804, 1))

    return input_image


def load_image_val(datapoint):
    input_image = tf.reshape(datapoint['ccf'], (804, 1))
    input_mask = tf.reshape(datapoint['rv'], (804, 1))

    return input_image, input_mask


# neural network functions
# 1D Convolution Block
def Conv1DBlock(input_tensor, n_filters, kernel_size=3, batchnorm=True):
    """
    add two 1D convolution layers to model

    Parameters
    ----------
    input_tensor: tensorflow compatible tensor for model fitting
    n_filters: output space dimensionality
    kernel_size: length of the 1D convolution window
    batchnorm: normalization layer (boolean)

    Returns
    -------

    """

    # first layer
    x = Conv1D(filters=n_filters, kernel_size=kernel_size, kernel_initializer='he_normal',
               padding='same')(input_tensor)
    if batchnorm:
        x = BatchNormalization()(x)
    x = Activation('relu')(x)

    # second layer
    x = Conv1D(filters=n_filters, kernel_size=kernel_size, kernel_initializer='he_normal',
               padding='same')(x)
    if batchnorm:
        x = BatchNormalization()(x)
    x = Activation('relu')(x)

    return x


# 1D Transposed Convolution
def Conv1DTranspose(input_tensor, filters, kernel_size, strides=2, padding='same'):
    """
    function to create a 1D transposed convolution
    based on the functionality of Keras Conv2DTranspose
    and their implementation of Conv1D

    Parameters
    ----------
    input_tensor: input tensor with the shape (batch_size, time_steps, dims)
    filters: output dimension
    kernel_size: size of the convolution kernel
    strides: convolution step size
    padding: same, valid

    Returns
    -------

    """

    x = Lambda(lambda x: BK.expand_dims(x, axis=2))(input_tensor)
    x = Conv2DTranspose(filters=filters, kernel_size=(kernel_size, 1), strides=(strides, 1), padding=padding)(x)
    x = Lambda(lambda x: BK.squeeze(x, axis=2))(x)
    return x


# U-NET Base Model
def get_unet(input_tensor, n_filters=16, dropout=0.1, batchnorm=True):
    """
    UNET basic 1D model for dealing with residual CCFs
    - tbh no idea if this is gonna work but we shall see

    Parameters
    ----------
    input_tensor: tensorflow compatible tensor for model fitting
    n_filters: output space dimensionality
    dropout: fraction of input units to drop
    batchnorm: normalization layer (boolean)

    Returns
    -------

    """

    # Contracting Path
    c1 = Conv1DBlock(input_tensor, n_filters * 1, kernel_size=3, batchnorm=batchnorm)
    p1 = MaxPooling1D((3))(c1)
    p1 = Dropout(dropout)(p1)

    c2 = Conv1DBlock(p1, n_filters * 2, kernel_size=3, batchnorm=batchnorm)
    p2 = MaxPooling1D((2))(c2)
    p2 = Dropout(dropout)(p2)

    c3 = Conv1DBlock(p2, n_filters * 4, kernel_size=3, batchnorm=batchnorm)
    p3 = MaxPooling1D((2))(c3)
    p3 = Dropout(dropout)(p3)

    # c4 = Conv1DBlock(p3, n_filters * 8, kernel_size=3, batchnorm=batchnorm)
    # p4 = MaxPooling1D((2))(c4)
    # p4 = Dropout(dropout)(p4)

    c5 = Conv1DBlock(p3, n_filters=n_filters * 16, kernel_size=3, batchnorm=batchnorm)

    # Expansive Path
    # u6 = Conv1DTranspose(c5, n_filters * 8, kernel_size=3, strides=(2), padding='same')
    # u6 = concatenate([u6, c4])
    # u6 = Dropout(dropout)(u6)
    # c6 = Conv1DBlock(u6, n_filters * 8, kernel_size=3, batchnorm=batchnorm)

    u7 = Conv1DTranspose(c5, n_filters * 4, kernel_size=3, strides=(2), padding='same')
    u7 = concatenate([u7, c3])
    u7 = Dropout(dropout)(u7)
    c7 = Conv1DBlock(u7, n_filters * 4, kernel_size=3, batchnorm=batchnorm)

    u8 = Conv1DTranspose(c7, n_filters * 2, kernel_size=3, strides=(2), padding='same')
    u8 = concatenate([u8, c2])
    u8 = Dropout(dropout)(u8)
    c8 = Conv1DBlock(u8, n_filters * 2, kernel_size=3, batchnorm=batchnorm)

    u9 = Conv1DTranspose(c8, n_filters * 1, kernel_size=3, strides=(3), padding='same')
    u9 = concatenate([u9, c1])
    u9 = Dropout(dropout)(u9)
    c9 = Conv1DBlock(u9, n_filters * 1, kernel_size=3, batchnorm=batchnorm)

    # use linear activation because this is a regression problem
    outputs = Conv1D(1, (1), activation='linear')(c9)
    model = Model(inputs=[input_tensor], outputs=[outputs])
    return model


# read in pickle
residual_pickle = os.path.join(CsvDir.CCFS, 'residual_ccfs.csv')
fpickle = pd.read_pickle(residual_pickle)
dates = fpickle.dates.values
residuals = fpickle.residual_ccf.values
rv_residual = fpickle.rv_residual.values
rv_error = fpickle.rv_error.values

# Set up velocity loop
config = Config.config
velocity_loop = np.arange(config['velocity_min'], config['velocity_max'], config['velocity_step']) + config['qrv']

# plug into the neural network and see what happens lol
model_h5 = 'map_unet_model.h5'
BATCH_SIZE = 16
EPOCHS = 10
train, test = train_test_split(dates, test_size=0.33, random_state=0)
train_inds = np.isin(dates, train)
test_inds = np.isin(dates, test)

### create tensorflow dataset
# convert to dictionaries
train_residuals = list(np.array(residuals)[train_inds])
train_velocity = [np.array(velocity_loop)]*len(train_residuals)
ccf_vel = np.array((train_residuals, train_velocity))
train_dict = {key: value for key, value in zip(["ccf_vel"], ccf_vel)}
train_dict['rv'] = list(rv_residual[train_inds])

val_residuals = list(np.array(residuals)[test_inds])
val_velocity = [np.array(velocity_loop)]*len(val_residuals)
ccf_vel_val = np.array((val_residuals, val_velocity))
val_dict = {key: value for key, value in zip(["ccf"], ccf_vel_val)}
val_dict['rv'] = list(rv_residual[test_inds])

# create tensorflow datasets
train_dataset = tf.data.Dataset.from_tensor_slices(train_dict)
val_dataset = tf.data.Dataset.from_tensor_slices(val_dict)


# reshape tensors to be 2D
# train = train_dataset.map(load_image_train, num_parallel_calls=tf.data.experimental.AUTOTUNE)
# val = val_dataset.map(load_image_val)

# train_dataset = train_dataset.prefetch(buffer_size=tf.data.experimental.AUTOTUNE)
# val_dataset = val.batch(BATCH_SIZE)

input_tensor = Input((804, 2), name='ccf_vel')
model = get_unet(input_tensor, n_filters=16, dropout=0.05, batchnorm=True)
model.compile(optimizer=Adam(), loss="binary_crossentropy", metrics=["accuracy"])

model.summary()

# callbacks
callbacks = [
    EarlyStopping(patience=10, verbose=1),
    ReduceLROnPlateau(factor=0.1, patience=5, min_lr=0.00001, verbose=1),
    ModelCheckpoint(model_h5, verbose=1, save_best_only=True, save_weights_only=True)
]

# fit model
results = model.fit(train_dataset, epochs=EPOCHS, callbacks=callbacks, validation_data=val_dataset)