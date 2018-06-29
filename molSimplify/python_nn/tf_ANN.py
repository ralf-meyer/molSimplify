# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## This script contains a neural network  ##########
#####  trained on octahedral metal-ligand          #######
########   bond distances and spin propensity  ###########
##########################################################


## import 
import keras
from keras.models import model_from_json
import numpy as np
import csv
from pkg_resources import resource_filename, Requirement
from molSimplify.Classes.globalvars import *
from molSimplify.python_nn.ANN import matrix_loader
import sys,os


## Functions
def data_rescale(scaled_dat,train_mean,train_var):
    d = np.shape(train_mean)[0]
    #print('unnormalizing with number of dimensions = ' +str(d))
    dat = (np.multiply(scaled_dat.T,np.sqrt(train_var),)+train_mean).T
    return(dat)
def data_normalize(data,train_mean,train_var):
    d = np.shape(train_mean)[0]
    #print('normalizing with number of dimensions = ' +str(d))
    scaled_dat = np.divide((data.T - train_mean),np.sqrt(train_var),).T
    return(scaled_dat)
    
def load_normalization_data(name):
    
    train_mean_x = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" + '/rescaling_data/'+name+'_mean_x.csv')
    with open(path_to_file,'r') as f:
        for lines in f.readlines():
            train_mean_x.append([float(lines.strip().strip('[]'))])
            
    train_var_x = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" + '/rescaling_data/'+name+'_var_x.csv')
    with open(path_to_file,'r') as f:
        for lines in f.readlines():
            train_var_x.append([float(lines.strip().strip('[]'))])

    train_mean_y = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" + '/rescaling_data/'+name+'_mean_y.csv')
    with open(path_to_file,'r') as f:
        for lines in f.readlines():
            train_mean_y.append([float(lines.strip().strip('[]'))])
    train_var_y = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" + '/rescaling_data/'+name+'_var_y.csv')
    with open(path_to_file,'r') as f:
        for lines in f.readlines():
            train_var_y.append([float(lines.strip().strip('[]'))])
    train_mean_x = np.array(train_mean_x)
    train_var_x = np.array(train_var_x)
    train_mean_y = np.array(train_mean_y)
    train_var_y = np.array(train_var_y)
    
    return train_mean_x,train_mean_y,train_var_x,train_var_y
    
    
def load_ANN_variables(predictor):
    if predictor in ['ls_ii','hs_ii','ls_iii','hs_iii']:
        key = 'geos/'+predictor+ '_vars'
    else:
        key = predictor+ '/'+predictor+ '_vars'
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" +key +'.csv')
    names = []
    with open(path_to_file,'r') as f:
        for lines in f.readlines():
            names.append(lines.strip())
    return names 

def load_training_data(predictor):
    if predictor in ['ls_ii','hs_ii','ls_iii','hs_iii']:
        key = 'geos/'+predictor+ '_bl_x'
    elif predictor == "split":
        key = predictor+ '/'+predictor+ '_x_41_OHE'
    else:
        key = predictor+ '/'+predictor+ '_x_OHE'
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" +key +'.csv')
    with open(path_to_file, "r") as f:
        csv_lines = list(csv.reader(f))
        #row_names = [row[0] for row in csv_lines]
        mat = [row for row in csv_lines[1:]]
    return mat
    
def load_keras_ann(predictor):
    ## this function loads the ANN for property
    ## "predcitor" 
    # disable TF output text to reduce console spam
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
    if predictor in ['ls_ii','hs_ii','ls_iii','hs_iii']:
        key = 'geos/'+predictor+ '_model'
    else:
        key = predictor+ '/model'
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" +key + '.json')
    json_file = open(path_to_file, 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    # load weights into  model
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" +key+ '.h5')
    loaded_model.load_weights(path_to_file)
    # complile model
    loaded_model.compile(loss="mse",optimizer='adam',
              metrics=['mse', 'mae', 'mape'])
    
    print("Keras/tf model loaded for " + str(predictor) + " from disk")
    return(loaded_model)


def tf_ANN_excitation_prepare(predictor,descriptors,descriptor_names):
    ## this function reforms the provided list of descriptors and their
    ## names to match the expectations of the target ANN model.
    ## it does NOT perfrom standardization
    
    print('preparing features for ' + str(predictor) +  ', recieved ' + str(len(descriptors)) + ' descriptors') 
    
    ## get variable names
    target_names = load_ANN_variables(predictor)
    print('model requires ' +  str(len(target_names)) + ' descriptors, attempting match')
    excitation = []
    valid = True
    for var_name in target_names:
        
        try:
            excitation.append(descriptors[descriptor_names.index(var_name)])
        except:
            print('looking for  ' + str(var_name))
            print('Error! variable  ' + str(var_name) + ' not found!')
            valid = False
            break
    excitation = np.array(excitation)
    excitation = np.reshape(excitation, (1,len(target_names)) )
    return excitation

def ANN_supervisor(predictor,descriptors,descriptor_names):
    print('ANN activated for ' + str(predictor)) 
    
    ## form the excitation in the corrrect order/variables
    excitation = tf_ANN_excitation_prepare(predictor,descriptors,descriptor_names)


    print('excitation is ' + str(excitation.shape))
    print('fetching non-dimensionalization data... ')
    train_mean_x,train_mean_y,train_var_x,train_var_y = load_normalization_data(predictor)
    print('rescaling input excitation...')
    excitation = data_normalize(excitation,train_mean_x,train_var_x)

    ## fetch ANN
    loaded_model = load_keras_ann(predictor)
    
    print('calling ANN model...')
    result = data_rescale(loaded_model.predict(excitation),train_mean_y,train_var_y)
    
    return result
         
def find_true_min_eu_dist(predictor,descriptors,descriptor_names):
    # returns scaled euclidean distance to nearest trainning 
    # vector in desciptor space
    
    ## form the excitation in the corrrect order/variables
    excitation = tf_ANN_excitation_prepare(predictor,descriptors,descriptor_names)
    
    ## getting train matrix info
    mat = load_training_data(predictor)
    train_mat = np.array(mat,dtype='float64')
    ## loop over rows
    min_dist = 1000
    min_ind = 0
    for i,rows in enumerate(train_mat):
        this_dist = np.linalg.norm(np.subtract(rows,np.array(excitation)))
        if this_dist < min_dist:
            min_dist = this_dist
            min_ind = i
            #best_row = rownames[i]
            min_row = rows
    
    # flatten min row
    min_row = np.reshape(min_row, excitation.shape) 
    
    #print('min dist is ' +str(min_dist) + ' at  ' + str(min_ind))
    
    # need to get normalized distances 
    train_mean_x,train_mean_y,train_var_x,train_var_y = load_normalization_data(predictor)
    scaled_excitation = data_normalize(excitation,train_mean_x,train_var_x)
    scaled_row = data_normalize(min_row,train_mean_x,train_var_x)
    min_dist = np.linalg.norm(np.subtract(scaled_row,(scaled_excitation)))
    
    return(min_dist)
    
    
    
    
