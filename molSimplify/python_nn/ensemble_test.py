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
from clf_analysis_tool import plot_scatter, plot_dist_err


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


def ensemble_maker(predictor, num=10):
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    mat = load_training_data(predictor)
    mat = np.array(mat, dtype='float64')
    train_mat = data_normalize(mat, train_mean_x, train_var_x)
    labels = load_training_labels(predictor)
    if not 'clf' in predictor:
        labels = np.array(labels, dtype='float64')
        labels = data_normalize(labels, train_mean_y, train_var_y)
    train_mat, labels = shuffle(train_mat, labels)
    train_mat = np.array_split(train_mat, num, axis=0)
    labels = np.array_split(labels, num, axis=0)
    info_dict = load_train_info(predictor)
    earlystop = EarlyStopping(monitor=info_dict['monitor'], min_delta=info_dict['min_delta'],
                              patience=info_dict['patience'],
                              verbose=1,
                              mode='auto')
    callbacks_list = [earlystop]
    for ii in range(num):
        train_feature = array_stack(train_mat, ii)
        train_labels = array_stack(labels, ii)
        loaded_model = load_keras_ann(predictor)
        current_model = reset_weights(loaded_model)
        history = current_model.fit(train_feature, train_labels,
                                    epochs=info_dict['epochs'], verbose=1,
                                    batch_size=info_dict['batch_size'],
                                    callbacks=callbacks_list)
        save_model(current_model, predictor, ii)


def ensemble_uq(predictor, descriptors=False, descriptor_names=False, suffix=False):
    key = get_key(predictor, suffix)
    base_path = resource_filename(Requirement.parse("molSimplify"), "molSimplify/tf_nn/" + key)
    base_path = base_path + 'ensemble_models'
    if not os.path.exists(base_path):
        print('Ensemble models do not exist now, training...')
        ensemble_maker(predictor)
    print('ANN activated for ' + str(predictor))
    model_list = glob.glob(base_path + '/*.h5')

    labels = load_training_labels(predictor)
    if not 'clf' in predictor:
        labels = np.array(labels, dtype='float64')
    else:
        labels = np.array(labels, dtype='int')
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    if (descriptors and descriptor_names):
        excitation = tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names)
        excitation = data_normalize(excitation, train_mean_x, train_var_x)
    else:
        mat = load_training_data(predictor)
        mat = np.array(mat, dtype='float64')
        train_mat = data_normalize(mat, train_mean_x, train_var_x)
        excitation = np.array(train_mat)
    print('excitation is ' + str(excitation.shape))
    print('actual label:', labels[:3])
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
        if not 'clf' in predictor:
            result = data_rescale(loaded_model.predict(excitation), train_mean_y, train_var_y)
        else:
            result = loaded_model.predict(excitation)
        results_list.append(result)
    results_list = np.transpose(np.array(results_list))
    result_mean, result_std = np.mean(results_list, axis=1), np.std(results_list, axis=1)
    labels = np.squeeze(labels, axis=1)
    error_for_mean = np.abs(labels - result_mean)
    return result_mean, result_std, error_for_mean


def mc_dropout_uq(predictor, descriptors=False, descriptor_names=False, num=100):
    labels = load_training_labels(predictor)
    if not 'clf' in predictor:
        labels = np.array(labels, dtype='float64')
    else:
        labels = np.array(labels, dtype='int')
    train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor)
    if (descriptors and descriptor_names):
        excitation = tf_ANN_excitation_prepare(predictor, descriptors, descriptor_names)
        excitation = data_normalize(excitation, train_mean_x, train_var_x)
    else:
        mat = load_training_data(predictor)
        mat = np.array(mat, dtype='float64')
        train_mat = data_normalize(mat, train_mean_x, train_var_x)
        excitation = np.array(train_mat)
    print('excitation is ' + str(excitation.shape))
    loaded_model = load_keras_ann(predictor)
    get_outputs = K.function([loaded_model.layers[0].input, K.learning_phase()],
                             [loaded_model.layers[-1].output])
    print('LOADED MODEL HAS ' + str(
        len(loaded_model.layers)) + ' layers, so latent space measure will be from first ' + str(
        len(loaded_model.layers) - 1) + ' layers')
    results_list = []
    err_list = []
    for ii in range(num):
        if not np.mod(ii, int(num / 10)):
            print('%d / %d' % (ii, num))
        if not 'clf' in predictor:
            results = data_rescale(np.array(get_outputs([excitation, 1])), train_mean_y,
                                   train_var_y)[0]
        else:
            results = np.array(get_outputs([excitation, 1]))[0]
        # print(results[:10], results.shape)
        # print(labels[:10], labels.shape)
        results = results.squeeze(axis=1)
        err = np.linalg.norm(labels - results) ** 2
        results_list.append(results)
        err_list.append(err)
    results_list = np.transpose(np.array(results_list))
    f = lambda tau: mc_dropout_logp(tau, np.array(err_list))
    tau = sp.optimize.minimize(f, 10).x
    result_mean, result_std = np.mean(results_list, axis=1), np.std(results_list, axis=1)
    result_std = np.sqrt(1 / tau + result_std ** 2)
    # print(tau, result_std[:3])
    labels = np.squeeze(labels, axis=1)
    error_for_mean = np.abs(labels - result_mean)
    return result_mean, result_std, error_for_mean


###########
predictor = 'geo_static_clf'
ensemble_maker(predictor, num=10)
# ensemble_uq(predictor)
# result_mean, result_std, error_for_mean = mc_dropout_uq(predictor)


# stds = [0.0075, 0.0125, 0.01875, 0.025, 0.0375, 0.05, 0.0625, 0.075, 0.0875, 0.1, 0.125, 0.15, 0.175, 0.2,
#         0.225, 0.25, 0.3, 0.375]
# plot_scatter(result_std, error_for_mean, xlabel='mc-dropout std', ylabel='prediction error',
#              figname='mc_dropout_var_err.pdf')
# plot_dist_err(pred_std=result_std, pred_err=error_for_mean, stds=stds,
#               label_x='mc-dropout std', label_y='accuracy',
#               lable_y2='ratio of data', figname='mc_dropout_var_acc_ratio.pdf')
