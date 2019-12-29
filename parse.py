import csv
import sys, getopt
import tensorflow as tf
from tensorflow import keras
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging


inputfile = sys.argv[2]
file = pd.read_csv(inputfile,header=None)
first_column = file.columns[0]
file = file.drop([first_column], axis=1)

result = file.values
print(result[0][0])
print(len(result))

# predictions = np.array([])
# inputfile = sys.argv[2]
# with open(inputfile) as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
#     line_count = 0
#     for row in csv_reader:        
#         line_count += 1
#         np.append(predictions, row)
#     # print(predictions)
#     print(f'Processed {line_count} lines.')
#     print(len(predictions))


#normalize in [0,1]

inputfile = tf.keras.utils.normalize(result, axis=1)
# print(inputfile)

# #sequential type of model
# model = tf.keras.models.Sequential()
# model.add(tf.keras.layers.Flatten())

# #N1 level (Fully Connected)
# model.add(tf.keras.layers.Dense(64, activation = tf.nn.relu))  

# Model
							# tf.logging.set_verbosity(tf.logging.ERROR)
logger = tf.get_logger()
logger.setLevel(logging.ERROR)

regressor = tf.contrib.learn.DNNRegressor(feature_columns=feature_cols, 
activation_fn = tf.nn.relu, hidden_units=[200, 100, 50, 25, 12])

# Reset the index of training
training_set.reset_index(drop = True, inplace =True)
def input_fn(data_set, pred = False):
    
    if pred == False:
        
        feature_cols = {k: tf.constant(data_set[k].values) for k in FEATURES}
        labels = tf.constant(data_set[LABEL].values)
        
        return feature_cols, labels

    if pred == True:
        feature_cols = {k: tf.constant(data_set[k].values) for k in FEATURES}
        
        return feature_cols

# Deep Neural Network Regressor with the training set which contain the data split by train test split
regressor.fit(input_fn=lambda: input_fn(training_set), steps=2000)

# Evaluation on the test set created by train_test_split
ev = regressor.evaluate(input_fn=lambda: input_fn(testing_set), steps=1)

# Display the score on the testing set
loss_score1 = ev["loss"]
print("Final Loss on the testing set: {0:f}".format(loss_score1))