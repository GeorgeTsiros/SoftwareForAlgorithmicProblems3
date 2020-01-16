import keras
import numpy as np
import pandas as pd
import sys
from keras import layers, optimizers, losses, metrics
from keras.models import load_model
import csv

# General Configurations
keras.backend.set_epsilon(1e-2)
np.set_printoptions(threshold=sys.maxsize)

# Read files
actualData = pd.read_csv('Datasets/actual.csv',header=None)

testFile = sys.argv[2]
testingData = pd.read_csv(testFile,header=None)

actualDataLabels = actualData.iloc[:, :1]
actualDataValues = actualData.iloc[:, 1:]
actualDataValues = actualDataValues.to_numpy()

testingDataLabels = testingData.iloc[:, :1]
testingDataValues = testingData.iloc[:, 1:]

model= load_model('Datasets/WindDenseNN.h5')
result = model.predict(testingDataValues)

# Generate MAE Loss
mae = keras.losses.MeanAbsoluteError()
lossmae = mae(actualDataValues, result)
print('MAE Loss: ', lossmae.numpy())  # Loss: 0.75

# Generate MSE Loss
mse = keras.losses.MeanSquaredError()
lossmse = mse(actualDataValues, result)
print('MSE Loss: ', lossmse.numpy())  # Loss: 0.75

# Generate MAPE Loss
mape = keras.losses.MeanAbsolutePercentageError()
lossmape = mape(actualDataValues, result)
print('MAPE Loss: ', lossmape.numpy())  # Loss: 5e+08

header = ['MAE: ' + str(lossmae.numpy()), ' MAPE: ' + str(lossmape.numpy()), ' MSE: ' + str(lossmse.numpy())]

combinedResults = np.concatenate((actualDataLabels, result), axis=1)
with open('predicted.csv', 'wt') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(header) # write header

pd.DataFrame(combinedResults).to_csv("predicted.csv", mode='a')