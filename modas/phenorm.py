import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import MinMaxScaler


def abundance_filter(d, abundance):
    d = d.loc[:,(d.var()-0) >= 1e-6]
    return d.loc[:, d.mean() >= abundance]


def missing_filter(d, missing_ratio):
    if d.applymap(np.isreal).all().sum() == d.shape[1]:
        d = d.loc[:, d.applymap(np.isnan).sum() <= d.shape[0] * missing_ratio]
        d = d.fillna(0)
    else:
        d = d.loc[:, d.applymap(np.isreal).sum() <= d.shape[0] * missing_ratio]
        d_array = d.values
        d_array[d.applymap(np.isreal)] = 0
        d.loc[:,:] = d_array
    d = d.loc[:, (d == 0).sum() <= d.shape[0] * missing_ratio]
    return d


def log2_scale(d):
    return np.log2(d+1)


def ln_scale(d):
    return np.log(d+1)


def log10_scale(d):
    return np.log10(d+1)


def normalize_scale(d):
    d = d + 1
    d = d.apply(lambda x: stats.boxcox(x)[0])
    d.loc[:,:] = MinMaxScaler().fit_transform(d.values)
    return d
