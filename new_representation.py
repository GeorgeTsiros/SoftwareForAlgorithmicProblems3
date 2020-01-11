import keras
import numpy as np
import pandas as pd
import sys
from keras import layers, optimizers, losses, metrics
from keras.models import load_model

# General Configurations
keras.backend.set_epsilon(1e-2)
np.set_printoptions(threshold=sys.maxsize, precision=16)

# Read files
actualData = pd.read_csv('Datasets/actual.csv')
testingData = pd.read_csv('Datasets/nn_representations.csv')

testingDataLabels = testingData.iloc[:, :1]
testingDataValues = testingData.iloc[:, 1:]

model= load_model('Datasets/WindDenseNN.h5')
model.summary()
weights = model.layers[0].get_weights()

newModel = keras.Sequential()
# Adds a densely-connected layer with 64 nodes to the model:
# input shape defines the number of input features (dimensions)
# activation defines the activation function of each layer
newModel.add(layers.Dense(64, activation='relu', input_shape=(128,)))
newModel.layers[0].set_weights(weights)
result = newModel.predict(testingDataValues)
np.savetxt('new_representation.csv', result, delimiter=',', fmt='%1.16f')   # X is an array