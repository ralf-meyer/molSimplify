# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## This script contains a neural network  ##########
#####  trained on octahedral metal-ligand          #######
########   bond distances and spin propensity  ###########
##########################################################


## import 
import keras
from keras import backend as K
from keras.models import model_from_json
from keras.optimizers import Adam
import numpy as np
import csv
from pkg_resources import resource_filename, Requirement
from molSimplify.Classes.globalvars import *
from molSimplify.python_nn.ANN import matrix_loader
import sys, os
import json


## Functions
def data_rescale(scaled_dat, train_mean, train_var):
    d = np.shape(train_mean)[0]
    # print('unnormalizing with number of dimensions = ' +str(d))
    dat = (np.multiply(scaled_dat.T, np.sqrt(train_var), ) + train_mean).T
    return (dat)


def data_normalize(data, train_mean, train_var):
    data = data.astype(float)  # Make sure the data is always in float form
    d = np.shape(train_mean)[0]
    # print('normalizing with number of dimensions = ' +str(d))
    scaled_dat = np.divide((data.T - train_mean), np.sqrt(train_var), ).T
    return (scaled_dat)


def load_normalization_data(name):
    train_mean_x = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),
                                     "molSimplify/tf_nn/" + '/rescaling_data/' + name + '_mean_x.csv')
    with open(path_to_file, 'r') as f:
        for lines in f.readlines():
            train_mean_x.append([float(lines.strip().strip('[]'))])

    train_var_x = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),
                                     "molSimplify/tf_nn/" + '/rescaling_data/' + name + '_var_x.csv')
    with open(path_to_file, 'r') as f:
        for lines in f.readlines():
            train_var_x.append([float(lines.strip().strip('[]'))])

    train_mean_y = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),
                                     "molSimplify/tf_nn/" + '/rescaling_data/' + name + '_mean_y.csv')
    with open(path_to_file, 'r') as f:
        for lines in f.readlines():
            train_mean_y.append([float(lines.strip().strip('[]'))])
    train_var_y = list()
    path_to_file = resource_filename(Requirement.parse("molSimplify"),
                                     "molSimplify/tf_nn/" + '/rescaling_data/' + name + '_var_y.csv')
    with open(path_to_file, 'r') as f:
        for lines in f.readlines():
            train_var_y.append([float(lines.strip().strip('[]'))])

    train_mean_x = np.array(train_mean_x)
    train_var_x = np.array(train_var_x)
    train_mean_y = np.array(train_mean_y)
    train_var_y = np.array(train_var_y)

    return train_mean_x, train_mean_y, train_var_x, train_var_y


def load_ANN_variables(predictor):
    if predictor in ['ls_ii', 'hs_ii', 'ls_iii', 'hs_iii']:
        key = 'geos/' + predictor + '_vars'
    elif predictor in ['homo', 'gap']:
        key = 'homolumo/' + predictor + '_vars'
    elif predictor in ['oxo', 'hat']:
        key = 'oxocatalysis/' + predictor + '_vars'
    else:
        key = predictor + '/' + predictor + '_vars'
    path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.csv')
    names = []
    with open(path_to_file, 'r') as f:
        for lines in f.readlines():
            names.append(lines.strip())
    return names


def load_training_data(predictor):
    if predictor in ['ls_ii', 'hs_ii', 'ls_iii', 'hs_iii']:
        key = 'geos/' + predictor + '_bl_x'
    elif predictor in ['homo', 'gap']:
        key = 'homolumo/' + predictor + '_train_x'
    elif predictor in ['oxo', 'hat']:
        key = 'oxocatalysis/' + predictor + '_train_x'
    elif predictor == "split":
        key = predictor + '/' + predictor + '_x_41_OHE'
    else:
        key = predictor + '/' + predictor + '_x_OHE'
    path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.csv')
    with open(path_to_file, "r") as f:
        csv_lines = list(csv.reader(f))
        # row_names = [row[0] for row in csv_lines]
        mat = [row for row in csv_lines[1:]]
    return mat

def load_test_data(predictor):
    if predictor in ['ls_ii','hs_ii','ls_iii','hs_iii']:
        key = 'geos/'+predictor+ '_bl_x' #Note, this test data is not available, will return train.
    elif predictor in ['homo','gap']:
        key = 'homolumo/'+predictor+'_test_x'
    elif predictor in ['oxo','hat']:
        key = 'oxocatalysis/'+predictor+ '_test_x'
    elif predictor == "split":
        key = predictor+ '/'+predictor+ '_x_41_OHE' #Note, this test data is not available, will return train
    else:
        key = predictor+ '/'+predictor+ '_x_OHE'
    path_to_file = resource_filename(Requirement.parse("molSimplify"),"molSimplify/tf_nn/" +key +'.csv')
    with open(path_to_file, "r") as f:
        csv_lines = list(csv.reader(f))
        #row_names = [row[0] for row in csv_lines]
        mat = [row for row in csv_lines[1:]]
    return mat

def load_training_labels(predictor):
    if predictor in ['ls_ii', 'hs_ii', 'ls_iii', 'hs_iii']:
        key = 'geos/' + predictor + '_bl_y'
    elif predictor in ['homo', 'gap']:
        key = 'homolumo/' + predictor + '_train_y'
    elif predictor in ['oxo', 'hat']:
        key = 'oxocatalysis/' + predictor + '_train_y'
    elif predictor == "split":
        key = predictor + '/' + predictor + '_y_41_OHE'
    else:
        key = predictor + '/' + predictor + '_y_OHE'
    path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.csv')
    with open(path_to_file, "r") as f:
        csv_lines = list(csv.reader(f))
        # row_names = [row[0] for row in csv_lines]
        mat = [row for row in csv_lines[1:]]
    return mat


def load_train_info(predictor):
    if predictor in ['ls_ii', 'hs_ii', 'ls_iii', 'hs_iii']:
        key = 'geos/' + predictor + '_info'
    elif predictor in ['homo', 'gap']:
        key = 'homolumo/' + predictor + '_info'
    elif predictor in ['oxo', 'hat']:
        key = 'oxocatalysis/' + predictor + '_info'
    else:
        print('!!!!! checkout!!!!!')
        key = predictor + '/model'
    path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.json')
    json_file = open(path_to_file, 'r')
    loaded_info_dict = json.loads(json_file.read())
    json_file.close()
    return loaded_info_dict


def load_keras_ann(predictor):
    ## this function loads the ANN for property
    ## "predcitor" 
    # disable TF output text to reduce console spam
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    if predictor in ['ls_ii', 'hs_ii', 'ls_iii', 'hs_iii']:
        key = 'geos/' + predictor + '_model'
    elif predictor in ['homo', 'gap']:
        key = 'homolumo/' + predictor + '_model'
    elif predictor in ['oxo', 'hat']:
        key = 'oxocatalysis/' + predictor + '_model'
    else:
        key = predictor + '/model'
    path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.json')
    json_file = open(path_to_file, 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    # load weights into  model
    path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.h5')
    loaded_model.load_weights(path_to_file)
    # complile model
    if predictor == 'homo':
        loaded_model.compile(loss="mse", optimizer=Adam(beta_2=1 - 0.0016204733101599046, beta_1=0.8718839135783554,
                                                        decay=7.770243145972892e-05, lr=0.0004961686075897741),
                             metrics=['mse', 'mae', 'mape'])
    elif predictor == 'gap':
        loaded_model.compile(loss="mse", optimizer=Adam(beta_2=1 - 0.00010929248596488832, beta_1=0.8406735969305784,
                                                        decay=0.00011224350434148253, lr=0.0006759924688701965),
                             metrics=['mse', 'mae', 'mape'])
    elif predictor in ['oxo', 'hat']:
        loaded_model.compile(loss="mse", optimizer=Adam(beta_2=0.9637165412871632, beta_1=0.7560951483268549,
                                                        decay=0.0006651401379502965, lr=0.0007727366541920176),
                             metrics=['mse', 'mae', 'mape'])
    else:
        loaded_model.compile(loss="mse", optimizer='adam',
                             metrics=['mse', 'mae', 'mape'])

    print("Keras/tf model loaded for " + str(predictor) + " from disk")

    return (loaded_model)


def tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names):
    ## this function reforms the provided list of descriptors and their
    ## names to match the expectations of the target ANN model.
    ## it does NOT perfrom standardization

    print('preparing features for ' + str(predictor) + ', recieved ' + str(len(descriptors)) + ' descriptors')

    ## get variable names
    target_names = load_ANN_variables(predictor)
    print('model requires ' + str(len(target_names)) + ' descriptors, attempting match')
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
    excitation = np.reshape(excitation, (1, len(target_names)))
    return excitation


def ANN_supervisor(predictor, descriptors, descriptor_names):
    print('ANN activated for ' + str(predictor))

    ## form the excitation in the corrrect order/variables
    excitation = tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names)

    print('excitation is ' + str(excitation.shape))
    print('fetching non-dimensionalization data... ')
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    print('rescaling input excitation...')

    excitation = data_normalize(excitation, train_mean_x, train_var_x)

    ## fetch ANN
    loaded_model = load_keras_ann(predictor)
    print('LOADED MODEL HAS ' + str(
        len(loaded_model.layers)) + ' layers, so latent space measure will be from first ' + str(
        len(loaded_model.layers) - 1) + ' layers')
    get_outputs = K.function([loaded_model.layers[0].input, K.learning_phase()],
                             [loaded_model.layers[len(loaded_model.layers) - 2].output])
    latent_space_vector = get_outputs([excitation, 0])  # Using test phase.

    print('calling ANN model...')
    result = data_rescale(loaded_model.predict(excitation), train_mean_y, train_var_y)
    return result, latent_space_vector


def find_true_min_eu_dist(predictor, descriptors, descriptor_names):
    # returns scaled euclidean distance to nearest trainning 
    # vector in desciptor space
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)

    ## form the excitation in the corrrect order/variables
    excitation = tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names)
    excitation = excitation.astype(float)  # ensure that the excitation is a float, and not strings
    scaled_excitation = data_normalize(excitation, train_mean_x, train_var_x)  # normalize the excitation
    ## getting train matrix info
    mat = load_training_data(predictor)
    train_mat = np.array(mat, dtype='float64')
    ## loop over rows
    min_dist = 100000000
    min_ind = 0
    for i, rows in enumerate(train_mat):
        scaled_row = np.squeeze(
            data_normalize(rows, train_mean_x.T, train_var_x.T))  # Normalizing the row before finding the distance
        this_dist = np.linalg.norm(np.subtract(scaled_row, np.array(scaled_excitation)))
        if this_dist < min_dist:
            min_dist = this_dist
            min_ind = i
            # best_row = rownames[i]
            min_row = rows

    # flatten min row
    min_row = np.reshape(min_row, excitation.shape)
    print('min dist is ' + str(min_dist) + ' at  ' + str(min_ind))
    if predictor in ['oxo', 'hat', 'homo', 'gap']:
        if predictor in ['homo', 'gap']:
            key = 'homolumo/' + predictor + '_train_names'
        elif predictor in ['oxo', 'hat']:
            key = 'oxocatalysis/' + predictor + '_train_names'
        path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.csv')
        with open(path_to_file, "r") as f:
            csv_lines = list(csv.reader(f))
            print('Closest Structure: ', csv_lines[min_ind + 1])
    # need to get normalized distances 

    ########################################################################################
    # Changed by Aditya on 08/13/2018. Previously, nearest neighbor was being found in the #
    # unnormalized space, and then that was normalized. This was resulting in bad nearest  #
    # neighbor candidate structures. Now routine normalizes before finding the distance.   #
    ########################################################################################

    # train_mean_x,train_mean_y,train_var_x,train_var_y = load_normalization_data(predictor)

    # scaled_excitation = data_normalize(excitation,train_mean_x,train_var_x)
    # scaled_row = data_normalize(min_row,train_mean_x,train_var_x)
    # min_dist = np.linalg.norm(np.subtract(scaled_row,(scaled_excitation)))
    return (min_dist)


def find_ANN_latent_dist(predictor, latent_space_vector):
    # returns scaled euclidean distance to nearest trainning 
    # vector in desciptor space
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)

    ## getting train matrix info
    mat = load_training_data(predictor)
    train_mat = np.array(mat, dtype='float64')
    ## loop over rows
    min_dist = 100000000
    min_ind = 0

    loaded_model = load_keras_ann(predictor)
    print('MEASURING LATENT SPACE DISTANCE!')
    print('LOADED MODEL HAS ' + str(
        len(loaded_model.layers)) + ' layers, so latent space measure will be from first ' + str(
        len(loaded_model.layers) - 1) + ' layers')
    get_outputs = K.function([loaded_model.layers[0].input, K.learning_phase()],
                             [loaded_model.layers[len(loaded_model.layers) - 2].output])

    for i, rows in enumerate(train_mat):
        scaled_row = np.squeeze(
            data_normalize(rows, train_mean_x.T, train_var_x.T))  # Normalizing the row before finding the distance
        latent_train_row = get_outputs([np.array([scaled_row]), 0])
        this_dist = np.linalg.norm(np.subtract(np.squeeze(latent_train_row), np.squeeze(latent_space_vector)))
        # print(this_dist)
        if this_dist < min_dist:
            min_dist = this_dist
            min_ind = i
            # best_row = rownames[i]
            min_row = rows

    # flatten min row
    print('min dist is ' + str(min_dist) + ' at  ' + str(min_ind))
    if predictor in ['oxo', 'hat', 'homo', 'gap']:
        if predictor in ['homo', 'gap']:
            key = 'homolumo/' + predictor + '_train_names'
        elif predictor in ['oxo', 'hat']:
            key = 'oxocatalysis/' + predictor + '_train_names'
        path_to_file = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key + '.csv')
        with open(path_to_file, "r") as f:
            csv_lines = list(csv.reader(f))
            print('Closest Structure: ', csv_lines[min_ind + 1])
    # need to get normalized distances 

    ########################################################################################
    # Changed by Aditya on 08/13/2018. Previously, nearest neighbor was being found in the #
    # unnormalized space, and then that was normalized. This was resulting in bad nearest  #
    # neighbor candidate structures. Now routine normalizes before finding the distance.   #
    ########################################################################################

    # train_mean_x,train_mean_y,train_var_x,train_var_y = load_normalization_data(predictor)

    # scaled_excitation = data_normalize(excitation,train_mean_x,train_var_x)
    # scaled_row = data_normalize(min_row,train_mean_x,train_var_x)
    # min_dist = np.linalg.norm(np.subtract(scaled_row,(scaled_excitation)))
    return (min_dist)


def save_model(model, predictor, num=None):
    if predictor in ['ls_ii', 'hs_ii', 'ls_iii', 'hs_iii']:
        key = 'geos/'
    elif predictor in ['homo', 'gap']:
        key = 'homolumo/'
    elif predictor in ['oxo', 'hat']:
        key = 'oxocatalysis/'
    else:
        key = predictor
    base_path = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key)
    base_path = base_path + 'ensemble_models'
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    if not num == None:
        name = '%s/%s_%d' % (base_path, predictor, num)
    else:
        name = '%s/%s' % (base_path, predictor)
    # serialize model to JSON
    model_json = model.to_json()
    with open("%s.json" % name, "w") as json_file:
        json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights("%s.h5" % name)
    print("Saved model !%s! to disk" % name.split('/')[-1])
