# -*- coding: utf-8 -*-
"""
#   HLApipeline#
#   Authors:    Limin Jiang (lxj423@med.miami.edu)
#   Date:       2024-3-1
#   Version:    0.2.0
# Default values
"""

import pyreadr
import numpy as np
import keras
import keras.backend as K
from keras import layers
import pandas as pd
import random
import tensorflow as tf
import os
np.random.seed(42)

def padding(array, xx, yy):
    """
    :param array: numpy array
    :param xx: desired height
    :param yy: desirex width
    :return: padded array
    """

    h = array.shape[0]
    w = array.shape[1]

    a = (xx - h) // 2
    aa = xx - a - h

    b = (yy - w) // 2
    bb = yy - b - w

    return np.pad(array, pad_width=((a, aa), (b, bb)), mode='constant')



def DataGenerator_temp(csvFile,idx_l0,idx_l1,data_L00,data_L11,data_L01,data_L10,Sizeof,sampleNum): 
    data_cluster1 = []
    data_cluster2 = []
    data_between1 = []
    data_between2 = []
    NUmber0 = round(len(idx_l0)*0.7)
    NUmber1 = round(len(idx_l1)*0.7)
    for ii in range(sampleNum):      
        
        index_t0 = random.sample(list(data_L00.index), NUmber0)
        index_t1 =  random.sample(list(data_L11.index), NUmber1)   
        data_L00_temp = data_L00.loc[index_t0,:].iloc[:,index_t0].to_numpy() 
        data_L00_temp = padding(data_L00_temp, Sizeof, Sizeof) 
        data_cluster1.append(data_L00_temp)
        data_L11_temp = data_L11.loc[index_t1,:].iloc[:,index_t1].to_numpy()
        data_L11_temp = padding(data_L11_temp, Sizeof, Sizeof) 
        data_cluster2.append(data_L11_temp)
        data_L01_temp = data_L01.loc[index_t0,:].iloc[:,index_t1].to_numpy()
        data_L01_temp = padding(data_L01_temp, Sizeof, Sizeof) 
        data_between1.append(data_L01_temp)
        data_L10_temp = data_L10.loc[index_t1,:].iloc[:,index_t0].to_numpy()
        data_L10_temp = padding(data_L10_temp, Sizeof, Sizeof) 
        data_between2.append(data_L10_temp)
        
        
    for ii in range(2*sampleNum):             
        
        index_t0 = random.sample(list(csvFile.index), NUmber0)
        index_t1 =  random.sample(list(csvFile.index), NUmber1)
        
        if (not set(index_t0).issubset(set(idx_l0))) and (not set(index_t1).issubset(set(idx_l1))):         
            data_L00_temp = csvFile.loc[index_t0,:].iloc[:,index_t0].to_numpy() 
            data_L00_temp = padding(data_L00_temp, Sizeof, Sizeof) 
            data_cluster1.append(data_L11_temp)
            data_L11_temp = csvFile.loc[index_t1,:].iloc[:,index_t1].to_numpy()
            data_L11_temp = padding(data_L11_temp, Sizeof, Sizeof) 
            data_cluster2.append(data_L00_temp)
            data_L01_temp = csvFile.loc[index_t0,:].iloc[:,index_t1].to_numpy()
            data_L01_temp = padding(data_L01_temp, Sizeof, Sizeof) 
            data_between1.append(data_L01_temp)
            data_L10_temp = csvFile.loc[index_t1,:].iloc[:,index_t0].to_numpy()
            data_L10_temp = padding(data_L10_temp, Sizeof, Sizeof)             
            data_between2.append(data_L10_temp)
    
    data_cluster1 = np.array(data_cluster1)
    data_cluster2 = np.array(data_cluster2)
    data_between1 = np.array(data_between1)
    data_between2 = np.array(data_between2)
    labels = np.concatenate([np.ones(sampleNum),np.zeros(len(data_cluster1) - sampleNum)])
    
    return data_cluster1,data_cluster2,data_between1,data_between2,labels
    
def contrastive_loss(alpha):
    def loss(y_true, y_pred):       
        anchor1, anchor2, positive, negative = y_pred[:, 0], y_pred[:, 1], y_pred[:, 2], y_pred[:, 3]
        pos_distance = K.sqrt(K.sum(K.square(anchor1 - positive), axis=-1))  
        neg_distance = K.sqrt(K.sum(K.square(anchor2 - negative), axis=-1))
        between_distance = K.sqrt(K.sum(K.square(negative - positive), axis=-1))
        loss_value = ( 1/2*(1 - y_true) * (pos_distance+neg_distance+between_distance) ) + ( 1/2*y_true*max(0,alpha - (pos_distance+neg_distance+between_distance)))
        return K.mean(loss_value)

    return loss



# Define Hyperparameters
num_classes = 2
batch_size = 128
epochs = 128
embedding_dim = 32
alpha = 0.5
sampleNum =1000



directory_path = os.getcwd()
all_files = os.listdir(directory_path)
rda_files = [file for file in all_files if  ".rda" in file]


for rda_file in rda_files:
    path_temp = os.path.join(directory_path, rda_file)
    csvFile = pyreadr.read_r(path_temp)

    sample_label = csvFile['sample_label']
    #sample_label = sample_label.sample(frac=1).reset_index(drop=True)
    csvFile = csvFile['data_all']
    csvFile =pd.concat([sample_label, csvFile], axis=1)

    
    if (len(csvFile) >= 20):            
        idx_l0 = csvFile.index[csvFile['sample_label'] == 0 ].tolist()
        idx_l1 = csvFile.index[csvFile['sample_label'] == 1 ].tolist()
        csvFile = csvFile.drop('sample_label',axis = 1)
        
        data_L00 = csvFile.loc[idx_l0,csvFile.columns[idx_l0]].reset_index(drop = True)
        data_L11 = csvFile.loc[idx_l1,csvFile.columns[idx_l1]].reset_index(drop = True)    
        data_L01 = csvFile.loc[idx_l0,csvFile.columns[idx_l1]].reset_index(drop = True)
        data_L10 = csvFile.loc[idx_l1,csvFile.columns[idx_l0]].reset_index(drop = True)
        
        
        
        
        Sizeof = max(len(idx_l0),len(idx_l1))    
        input_shape_temp = (Sizeof, Sizeof, 1) 
        
        
        
        data_cluster1,data_cluster2,data_between1,data_between2,y_labels = DataGenerator_temp(csvFile,idx_l0,idx_l1,data_L00,data_L11,data_L01,data_L10,Sizeof,sampleNum)
        
        
        
        anchor_input1 = layers.Input(shape=input_shape_temp, name="anchor_input1")
        anchor_input2 = layers.Input(shape=input_shape_temp, name="anchor_input2")
        positive_input = layers.Input(shape=input_shape_temp, name="positive_input")
        negative_input = layers.Input(shape=input_shape_temp, name="negative_input")
        
        
        encoder = keras.Sequential(
            [
                layers.Conv2D(16, (4, 4), activation="relu", input_shape=input_shape_temp),
                layers.MaxPooling2D((2, 2)),
                layers.Conv2D(32, (4, 4), activation="relu"),
                layers.MaxPooling2D((2, 2)),
                layers.Flatten(),
                layers.Dense(embedding_dim, activation="sigmoid"),
                layers.LayerNormalization(axis=-1),
            ],
            name="encoder",
        )
        
        
        
        encoded_anchor1 = encoder(anchor_input1)
        encoded_anchor2 = encoder(anchor_input2)
        encoded_positive = encoder(positive_input)
        encoded_negative = encoder(negative_input)
        
        merged_output = layers.concatenate([encoded_anchor1, encoded_anchor2, encoded_positive, encoded_negative], axis=-1, name="merged_layer")
        
        model = keras.Model(inputs=[anchor_input1, anchor_input2, positive_input, negative_input], outputs=merged_output, name="triplet_model")
        
        
        class DataGenerator(keras.utils.Sequence):
            def __init__(self, x, y, z, a, y_labels,batch_size,num_classes, alpha):
                self.x = x
                self.y = y
                self.z = z
                self.a = a
                self.y_labels = y_labels
                self.batch_size = batch_size
                self.num_classes = num_classes
                self.alpha = alpha
        
            def __len__(self):
                return int(np.ceil(len(self.x)) / float(self.batch_size))
        
            def __getitem__(self, index):
                anchor1 = self.x[index * self.batch_size: (index + 1) * self.batch_size]
                anchor2 = self.x[index * self.batch_size: (index + 1) * self.batch_size]
                cluster1 = self.y[index * self.batch_size: (index + 1) * self.batch_size]
                cluster2 = self.z[index * self.batch_size: (index + 1) * self.batch_size]        
        
                return [anchor1,anchor2, cluster1, cluster2], np.zeros((self.batch_size,))
        
        generator = DataGenerator(data_between1.reshape(-1, Sizeof, Sizeof, 1), data_between2.reshape(-1, Sizeof, Sizeof, 1), data_cluster1.reshape(-1, Sizeof, Sizeof, 1), data_cluster2.reshape(-1, Sizeof, Sizeof, 1),y_labels, batch_size,num_classes, alpha)
        
        learning_rate = 0.0000001
        momentum = 0.5
        nesterov = True
        optimizer = tf.keras.optimizers.SGD(learning_rate=learning_rate, momentum=momentum, nesterov=nesterov)
        
        
        model.compile(loss=contrastive_loss(alpha), optimizer= optimizer,run_eagerly=True)
        model.fit(generator, epochs=epochs)
        
        
        
        
        encoder = model.get_layer("encoder")
        embeddings_cluster1_train = encoder.predict(data_cluster1.reshape(-1, Sizeof, Sizeof, 1))
        embeddings_cluster2_train = encoder.predict(data_cluster2.reshape(-1, Sizeof, Sizeof, 1))
        embeddings_between1_train = encoder.predict(data_between1.reshape(-1, Sizeof, Sizeof, 1))
        embeddings_between2_train = encoder.predict(data_between2.reshape(-1, Sizeof, Sizeof, 1))
        
        
        
        
        
        
        data_cluster1_test,data_cluster2_test,data_between1_test,data_between2_test,y_labels_test = DataGenerator_temp(csvFile,idx_l0,idx_l1,data_L00,data_L11,data_L01,data_L10,Sizeof,10)
        
        embeddings_cluster1_test = encoder.predict(data_cluster1_test.reshape(-1, Sizeof, Sizeof, 1))
        embeddings_cluster2_test = encoder.predict(data_cluster2_test.reshape(-1, Sizeof, Sizeof, 1))
        embeddings_between1_test = encoder.predict(data_between1_test.reshape(-1, Sizeof, Sizeof, 1))
        embeddings_between2_test = encoder.predict(data_between2_test.reshape(-1, Sizeof, Sizeof, 1))
        
        
        inputs_train=np.concatenate((embeddings_cluster1_train, embeddings_cluster2_train, embeddings_between1_train, embeddings_between2_train), axis=1)
        inputs_test=np.concatenate((embeddings_cluster1_test, embeddings_cluster2_test, embeddings_between1_test, embeddings_between2_test), axis=1)
        
        
        mlp_model = keras.Sequential([
            layers.Dense(64, activation="relu", input_shape=(embedding_dim*4,)),
            layers.Dense(1, activation="sigmoid")
        ], name="mlp_model")
        
        # 编译并训练MLP
        mlp_model.compile(loss="binary_crossentropy", optimizer="adam", metrics=["accuracy"])
        mlp_model.fit(inputs_train, y_labels, batch_size=batch_size, epochs=epochs, validation_data=(inputs_test, y_labels_test))
        
        test_loss, test_acc = mlp_model.evaluate(inputs_test, y_labels_test)
        print("Test accuracy:", test_acc)
        
        
        embeddings_data_L00 = encoder.predict(padding(data_L00, Sizeof, Sizeof).reshape(-1, Sizeof, Sizeof, 1))
        embeddings_data_L11  = encoder.predict(padding(data_L11, Sizeof, Sizeof).reshape(-1, Sizeof, Sizeof, 1))  
        embeddings_data_L01  = encoder.predict(padding(data_L01, Sizeof, Sizeof).reshape(-1, Sizeof, Sizeof, 1))
        embeddings_data_L10  = encoder.predict(padding(data_L10, Sizeof, Sizeof).reshape(-1, Sizeof, Sizeof, 1))
        
        
        Values = mlp_model.predict(np.concatenate((embeddings_data_L00,embeddings_data_L11,embeddings_data_L01,embeddings_data_L10), axis=1))
        
        
        
        my_string = path_temp + "   " +  str(Values) + "\n"
        
        file_path = r'Results.txt' 
        with open(file_path, 'a+') as file:
            file.write(my_string)
        
        
        
        
        df = pd.DataFrame({'L00': embeddings_data_L00.flatten(), 'L11': embeddings_data_L11.flatten(), 'L01': embeddings_data_L01.flatten(), 'L10': embeddings_data_L10.flatten()})
        file_path = path_temp.replace(".rda", ".csv")
        df.to_csv(file_path, index=False)
        
        
        
        
        
        
        
        
        
        
