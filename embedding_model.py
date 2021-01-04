# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:37:34 2020

@author: harsh
"""

import numpy as np
from tensorflow import keras
from keras import layers
import pandas as pd
from keras.utils.np_utils import to_categorical
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Flatten, Activation
from keras.optimizers import Adam, RMSprop
from sklearn.metrics import precision_recall_curve
from keras import initializers
from keras import backend as K
from numpy.random import seed
import tensorflow as tf
from sklearn.decomposition import PCA
from scipy import stats

expression = pd.read_csv("Gene2vec/pre_trained_expression.csv")
embeddings = pd.read_csv("Gene2vec/pre_trained_embeddings.csv", header = None)
embeddings = embeddings.iloc[:,2:202]
expression_data = expression.iloc[:,7:expression.shape[1]]
labels = expression.iloc[:,4]
labels = to_categorical(labels.factorize()[0])




k = 5
n = 32
genes = 11751

#PCA
pca = PCA(n_components = n)
components = pd.DataFrame(pca.fit_transform(embeddings))


#leave one out CV ATTENTION
index = -1    
pred_probs = np.empty((1,k))
attention_scores = np.empty((195, genes))

seed(123)
def my_init(shape, dtype=None):
    return K.variable(value = components, dtype = dtype)

for i in range(1,196):
    print(i)

    train_data = expression_data.drop(expression_data.index[index+1])
    test_data = expression_data.iloc[index+1:index+2,:]
    train_labels = np.delete(labels, index+1, axis = 0)
    model = tf.keras.Sequential()
    #model.add(tf.keras.layers.Dropout(0.2, input_shape = (genes,)))
    model.add(tf.keras.layers.Dense(32, input_dim = genes, use_bias = False))
    model.add(tf.keras.layers.Dense(genes, activation = "softmax", input_dim = 32, use_bias = True))
    model.add(tf.keras.layers.Dense(n, use_bias = True, kernel_initializer = my_init, bias_initializer = tf.keras.initializers.Zeros(), trainable = True))
    model.add(tf.keras.layers.Dense(k,activation = "softmax"))
    
    model.compile(loss = 'categorical_crossentropy', optimizer = tf.keras.optimizers.Adam(lr = 0.02), metrics = ['accuracy'])
    model.fit(train_data, train_labels, epochs = 75, batch_size = 10, shuffle = True)
    pred_probs = np.append(pred_probs, model.predict(test_data), axis = 0)
    extractor = keras.Model(inputs=model.input,
                        outputs=[model.layers[2].output])
    output = np.array(memoryview(tf.constant(extractor(test_data))))
    attention_scores[i-1] = output
    index = index + 1

    
np.savetxt("attention_predprobs.csv", pred_probs, delimiter = " ")

np.savetxt("attention_data.csv", attention_scores, delimiter = ",")


#leave one out CV NON-ATTENTION
index = -1    
pred_probs = np.empty((1,k))


seed(123)
def my_init(shape, dtype=None):
    return K.variable(value = embeddings, dtype = dtype)

for i in range(1,196):
    print(i)

    train_data = expression_data.drop(expression_data.index[index+1])
    test_data = expression_data.iloc[index+1:index+2,:]
    train_labels = np.delete(labels, index+1, axis = 0)
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(n, input_dim = 11751, activation = "relu",use_bias = True, kernel_initializer = my_init, bias_initializer = tf.keras.initializers.Zeros(), trainable = True))
    model.add(tf.keras.layers.Dense(k,activation = "softmax"))
    
    model.compile(loss = 'categorical_crossentropy', optimizer = keras.optimizers.Adam(lr = 0.00263), metrics = ['accuracy'])
    model.fit(train_data, train_labels, epochs = 75, batch_size = 10, shuffle = True)
    pred_probs = np.append(pred_probs, model.predict(test_data), axis = 0)
    
    index = index + 1
