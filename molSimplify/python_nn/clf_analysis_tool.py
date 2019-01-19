import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances, euclidean_distances
from keras import backend as K


def plot_scatter(x, y,
                 xlabel=False, ylabel=False,
                 show=False, figname=False):
    try:
        plt.style.use('myline')
    except:
        pass
    fig = plt.figure(figsize=(8, 6))
    fig.add_subplot(111)
    plt.scatter(x, y)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    plt.tight_layout()
    if show:
        plt.show()
    if figname:
        fig.savefig(figname)


def plot_scatter_colored(x_axis, y_axis, c_axis,
                         xlabel, ylabel,
                         legend=None,
                         figname='tmp.pdf'):
    try:
        plt.style.use('myline')
    except:
        pass
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(x_axis, y_axis, c=c_axis,
                alpha=0.5, cmap='cool')
    cb = plt.colorbar()
    cb.set_label('Error')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if not legend == None:
        plt.legend(legend)
    plt.tight_layout()
    plt.show()
    if figname:
        fig.savefig(figname)


def simp_plot_line(y_axis, x_axis=None,
                   label_x='x_axis', label_y='y_axis',
                   log_x=False, log_y=False,
                   figname='tmp.pdf',
                   show=False, line=False):
    try:
        plt.style.use('myline')
    except:
        pass
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    if not line:
        ax.plot(x_axis, y_axis, '*')
    else:
        ax.plot(x_axis, y_axis, '-')
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    plt.tight_layout()
    if show:
        plt.show()
    if figname:
        fig.savefig(figname)


def plot_double_line(x_axis, y1_axis, y2_axis,
                     label_x='x_axis', label_y='y_axis',
                     lable_y2='y_axis1',
                     log_x=False, log_y=False,
                     legend=None,
                     figname='tmp.pdf',
                     show=False, y1_err_flag=False, y1_err_bar=None):
    try:
        plt.style.use('myline')
    except:
        pass
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    if y1_err_flag:
        p1 = ax.errorbar(x_axis, y1_axis, yerr=y1_err_bar, fmt='-bs')
    else:
        p1 = ax.plot(x_axis, y1_axis, '-bs')
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax2 = ax.twinx()
    p2 = ax2.plot(x_axis, y2_axis, '--r*')
    ax2.set_ylabel(lable_y2, color='r')
    ax2.tick_params('y', colors='r')
    plt.tight_layout()
    if not legend == None:
        plt.legend((p1[0], p2[0]), legend)
    if show:
        plt.show()
    if figname:
        fig.savefig(figname)


def get_acc(pred_std, pred_err, stds):
    pred_err_arr = []
    for target_std in stds:
        _pred_err = []
        for idx, _std in enumerate(pred_std):
            if _std < target_std:
                _pred_err.append(pred_err[idx])
        pred_err_arr.append(_pred_err)
    acc, auc, ratio = [], [], []
    for idx, _ in enumerate(pred_err_arr):
        pred_err_now = pred_err_arr[idx]
        acc_arr = [1 if pred_err_now[ii] < 0.5 else 0 for ii in range(len(pred_err_now))]
        acc.append(np.mean(acc_arr))
        ratio.append(1.0 * len(pred_err_now) / len(pred_std))
    print('stds', stds)
    print('accuracy', acc)
    print('ratio', ratio)
    return stds, np.array(acc), np.array(ratio)


def plot_dist_err(pred_std, pred_err, stds, label_x=False, label_y=False,
                  lable_y2=False, figname=False):
    std_arr, acc_arr, ratio_arr = get_acc(pred_std, pred_err, stds)
    plot_double_line(x_axis=std_arr, y1_axis=acc_arr, y2_axis=ratio_arr,
                     label_x=label_x, label_y=label_y,
                     lable_y2=lable_y2,
                     log_x=False, log_y=False,
                     legend=None,
                     figname=figname,
                     show=True)


def plot_metrics_correlation(metric1, metric2, pred_err,
                             xlabel=False, ylabel=False, figname='dist_relations.pdf'):
    match = [0 if x < 0.5 else 1 for x in pred_err]
    plot_scatter_colored(x_axis=metric1, y_axis=metric2, c_axis=match,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         legend=None,
                         figname=figname)


def dist_neighbor(fmat1, fmat2, labels, l=5, dist_ref=1):
    dist_mat = pairwise_distances(fmat1, fmat2, 'manhattan')
    dist_mat = dist_mat * 1.0 / dist_ref
    dist_avrg, dist_list, labels_list = [], [], []
    # print('shape of dist_mat:', dist_mat.shape)
    for ele in dist_mat:
        dist_arr = np.round(np.array(ele), 4)
        if not dist_ref == 1:
            _count = (dist_arr < 10).sum()
            _count = l if _count < l else _count
            _count = _count if _count < 300 else 300
        else:
            _count = l
        ind = dist_arr.argsort()[:_count]
        _dist = dist_arr[ind]
        dist_list.append(_dist)
        _labels = np.array([labels[x][0] for x in ind])
        labels_list.append(_labels)
        if _dist.all() > 1e-4:
            dist_avrg.append(np.mean(_dist[:l]))
        else:
            dist_avrg.append(np.mean(_dist[:l]) * float(l) / (l - 1))
    # print('-----mean: %f, std: %f---' % (np.mean(dist_avrg), np.std(dist_avrg)))
    dist_avrg = np.array(dist_avrg)
    dist_list = np.array(dist_list)
    labels_list = np.array(labels_list)
    return dist_avrg, dist_list, labels_list


def array_stack(array, _idx):
    for ii in range(len(array)):
        if not ii == _idx:
            if 'out_arr' not in dir():
                out_arr = array[ii]
            else:
                out_arr = np.concatenate((out_arr, array[ii]))
    return out_arr


def dist_penalty(d):
    return np.exp(-1 * d ** 2)


def get_entropy(dists, neighbor_targets):
    entropies = []
    for ii, _neighbor_targets in enumerate(neighbor_targets):
        p0, p1 = dist_penalty(2), dist_penalty(2)
        for idx, tar in enumerate(_neighbor_targets):
            tar = int(tar)
            d = dists[ii][idx]
            if d <= 10:
                if d != 0:
                    if tar == 0:
                        p0 += dist_penalty(d)
                    elif tar == 1:
                        p1 += dist_penalty(d)
                else:
                    if tar == 0:
                        p0 += 100
                    elif tar == 1:
                        p1 += 100
            _sum = p0 + p1
        p0 = p0 / _sum
        p1 = p1 / _sum
        if p1 == 0 or p0 == 0:
            entropies.append(0)
        else:
            entropies.append(-(p0 * np.log(p0) + p1 * np.log(p1)))
    return np.array(entropies)


def get_layer_outputs(model, layer_index, input,
                      training_flag=False):
    get_outputs = K.function([model.layers[0].input, K.learning_phase()],
                             [model.layers[layer_index].output])
    nn_outputs = get_outputs([input, training_flag])[0]
    return nn_outputs


def lse_trust(lse):
    level = False
    if lse < 0.2:
        level = 'very high'
    elif lse < 0.3:
        level = 'high'
    elif lse < 0.4:
        level = 'medium'
    elif lse < 0.5:
        level = 'low'
    else:
        level = 'very low'
    return level
