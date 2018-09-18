from keras import backend as K
from keras.callbacks import EarlyStopping
from keras.models import model_from_json
from tf_ANN import *
from sklearn.utils import shuffle
from pkg_resources import resource_filename, Requirement
import glob
import numpy as np
import scipy as sp
import math


def reset_weights(model):
    session = K.get_session()
    for layer in model.layers:
        if hasattr(layer, 'kernel_initializer'):
            layer.kernel.initializer.run(session=session)
    return model


def array_stack(array, _idx):
    for ii in range(len(array)):
        if not ii == _idx:
            if 'out_arr' not in dir():
                out_arr = array[ii]
            else:
                out_arr = np.concatenate((out_arr, array[ii]))
    return out_arr


def mc_dropout_logp(tau, err):
    T = err.shape[0]
    logp = sp.misc.logsumexp(-0.5 * tau * err)
    logp -= np.log(T)
    logp -= 0.5 * np.log(np.power(tau, -1))
    logp -= 0.5 * np.log(2 * math.pi)
    return (logp)
    
def ensemble_maker_inner(train_mat,labels,model_gen_function, info_dict,num=10):
    ## contains core functions to make ensemble models
    ## from training data and labels
    ## model_gen_function is a functiont that takes NO arguments and returns a keras model
    ## info_dict is a dictionary of training info 
    train_mat, labels = shuffle(train_mat, labels)
    train_mat = np.array_split(train_mat, num, axis=0)
    labels = np.array_split(labels, num, axis=0)
    earlystop = EarlyStopping(monitor=info_dict['monitor'], min_delta=info_dict['min_delta'],
                              patience=info_dict['patience'],
                              verbose=0,
                              mode='auto')
    callbacks_list = [earlystop]
    model_list = []
    for ii in range(num):
        train_feature = array_stack(train_mat, ii)
        train_labels = array_stack(labels, ii)
        loaded_model = model_gen_function() # note the call to gen new model
        current_model = reset_weights(loaded_model)
        history = current_model.fit(train_feature, train_labels,
                                    epochs=info_dict['epochs'], verbose=0,
                                    batch_size=info_dict['batch_size'],
                                    callbacks=callbacks_list)
        model_list.append(current_model)
    return(model_list)

def ensemble_maker(predictor, num=10):
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    mat = load_training_data(predictor)
    mat = np.array(mat, dtype='float64')
    train_mat = data_normalize(mat, train_mean_x, train_var_x)
    labels = load_training_labels(predictor)
    labels = np.array(labels, dtype='float64')
    labels = data_normalize(labels, train_mean_y, train_var_y)
    info_dict = load_train_info(predictor)
    model_list = ensemble_maker_inner(train_mat=train_mat,
                                      labels=labels,
                                      model_gen_function = lambda : load_keras_ann(predictor),
                                      info_dict = info_dict,num=num)
    for ii,current_model in enumerate(model_list):
        save_model(current_model, predictor, ii)


def ensemble_uq(predictor, descriptors, descriptor_names):
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
        print('Ensemble models do not exist now, training...')
        ensemble_maker(predictor)
    print('ANN activated for ' + str(predictor))

    ## form the excitation in the corrrect order/variables
    # excitation = tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names)

    model_list = glob.glob(base_path + '/*.h5')

    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    excitation = tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names)
    excitation = data_normalize(excitation, train_mean_x, train_var_x)
    print('excitation is ' + str(excitation.shape))
    print('fetching non-dimensionalization data... ')
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    print('rescaling input excitation...')
    labels = load_training_labels(predictor)
    labels = np.array(labels, dtype='float64')
    print('actrual label:', labels[:3])

    excitation = data_normalize(excitation, train_mean_x, train_var_x)
    results_list = []
    for idx, model in enumerate(model_list):
        _base = model.strip('.h5')
        json_file = open(_base + '.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        # load weights into  model
        loaded_model.load_weights(model)
        # complile model
        loaded_model.compile(loss="mse", optimizer='adam',
                             metrics=['mse', 'mae', 'mape'])
        result = data_rescale(loaded_model.predict(excitation), train_mean_y, train_var_y)
        results_list.append(result)
    results_list = np.array(results_list)
    result_mean, result_std = np.mean(results_list, axis=0), np.std(results_list, axis=0)
    return result_mean, result_std


def mc_dropout_uq(predictor, descriptors, descriptor_names, num=100):
    excitation = tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names)
    labels = load_training_labels(predictor)
    labels = np.array(labels, dtype='float64')
    # train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    # mat = load_training_data(predictor)
    # mat = np.array(mat, dtype='float64')
    # train_mat = data_normalize(mat, train_mean_x, train_var_x)
    excitation = np.array(train_mat)
    print('excitation is ' + str(excitation.shape))
    print('fetching non-dimensionalization data... ')
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    print('rescaling input excitation...')

    excitation = data_normalize(excitation, train_mean_x, train_var_x)

    ## fetch ANN
    loaded_model = load_keras_ann(predictor)
    get_outputs = K.function([loaded_model.layers[0].input, K.learning_phase()],
                             [loaded_model.layers[-1].output])
    print('LOADED MODEL HAS ' + str(
        len(loaded_model.layers)) + ' layers, so latent space measure will be from first ' + str(
        len(loaded_model.layers) - 1) + ' layers')
    results_list = []
    err_list = []
    for ii in range(num):
        results = data_rescale(np.array(get_outputs([excitation, 1])), train_mean_y,
                               train_var_y)[0]
        results = results.squeeze(axis=1)
        err = np.linalg.norm(labels - results) ** 2
        results_list.append(results)
        err_list.append(err)
    results_list = np.array(results_list)
    f = lambda tau: mc_dropout_logp(tau, np.array(err_list))
    tau = sp.optimize.minimize(f, 10).x
    result_mean, result_std = np.mean(results_list, axis=0), np.std(results_list, axis=0)
    result_std = np.sqrt(1 / tau + result_std ** 2)
    print(result_std[:3])
    return result_mean, result_std


###########
predictor = 'homo'
# ensemble_maker(predictor, num=3)
ensemble_uq(predictor)
# mc_dropout_uq(predictor)
