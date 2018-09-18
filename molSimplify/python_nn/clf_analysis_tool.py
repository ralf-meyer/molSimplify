import matplotlib.pyplot as plt
import numpy as np


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
