import keras
import numpy as np
import pandas as pd
import sys
from keras import layers, optimizers, losses, metrics
from keras.models import load_model

# General Configurations
keras.backend.set_epsilon(1e-2)
np.set_printoptions(threshold=sys.maxsize)

# Read files
inputFile = sys.argv[2]
actualData = pd.read_csv(inputFile,header=None)

testFile = sys.argv[2]
testingData = pd.read_csv(testFile,header=None)

actualDataLabels = actualData.iloc[:, :1]
actualDataValues = actualData.iloc[:, 1:]
actualDataValues = actualDataValues.to_numpy()

testingDataLabels = testingData.iloc[:, :1]
testingDataValues = testingData.iloc[:, 1:]

result = pretrainedmodel.predict(testingDataValues)
# print to CSV instead of STDOUT
# print(result)

# Generate MAE Loss
mae = keras.losses.MeanAbsoluteError()
loss = mae(actualDataValues, result)
print('MAE Loss: ', loss.numpy())  # Loss: 0.75

# Generate MSE Loss
mse = keras.losses.MeanSquaredError()
loss = mse(actualDataValues, result)
print('MSE Loss: ', loss.numpy())  # Loss: 0.75

# Generate MAPE Loss
mape = keras.losses.MeanAbsolutePercentageError()
loss = mape(actualDataValues, result)
print('MAPE Loss: ', loss.numpy())  # Loss: 5e+08

# actualDataValuesColumns = list()
# for column in actualDataValues.T:
#     actualDataValuesColumns.append(column)

# resultColumns = list()
# for column in result.T:
#     resultColumns.append(column)