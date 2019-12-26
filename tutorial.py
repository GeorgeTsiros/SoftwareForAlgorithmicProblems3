# link to tutorial: https://www.youtube.com/watch?v=wQ8BIBpya2k&feature=emb_title

# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt

print(tf.__version__)


mnist = tf.keras.datasets.mnist
(x_train, y_train),(x_test, y_test) = mnist.load_data()
x_train = tf.keras.utils.normalize(x_train, axis=1)
x_test = tf.keras.utils.normalize(x_test, axis=1)

#sequential type of model
model = tf.keras.models.Sequential()

# input layer - flatten the image
model.add(tf.keras.layers.Flatten())

# 2 hidden layers with 128 neurons, activated with rectified linear function

# first layer (default) with activation function
model.add(tf.keras.layers.Dense(128, activation = tf.nn.relu))

# second layer 
model.add(tf.keras.layers.Dense(128, activation = tf.nn.relu))

#output layer with 10 classifications with probability distribution (softmax)
model.add(tf.keras.layers.Dense(10, activation = tf.nn.softmax))

#parameters of training
model.compile(optimizer = 'adam',loss = 'sparse_categorical_crossentropy', metrics = ['accuracy'])

#train the modell
model.fit(x_train, y_train, epochs = 3)

#evaluate the model
val_loss,val_accuracy = model.evaluate (x_test, y_test);
print(val_loss)
print(val_accuracy)

#save/load model
model.save('test.model')
new_model = tf.keras.models.load_model('test.model')

#make predictions
predictions = new_model.predict(x_test)
print(np.argmax(predictions[0]))




# print(x_train[0])
# plt.imshow(x_train[0])
# plt.show()