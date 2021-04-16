import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import norm
from sklearn.preprocessing import MinMaxScaler
import warnings

def abundance_filter(d, abundance):
    return d.loc[:, d.mean() >= abundance]


def isDigit(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def missing_filter(d, missing_ratio):
    if d.applymap(np.isreal).all().sum() == d.shape[1]:
        d = d.loc[:, d.applymap(np.isnan).sum() <= d.shape[0] * missing_ratio]
        d = d.fillna(0)
    else:
        d = d.loc[:, d.applymap(lambda x:isDigit(x)).sum() >= d.shape[0] * missing_ratio]
        d_array = d.values
        d_array[~d.applymap(lambda x:isDigit(x))] = 0
        d.loc[:, :] = d_array
        d = d.astype(float)
    #d = d.loc[:, (d == 0).sum() <= d.shape[0] * missing_ratio]
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


def ppoints(n, a=None):
    try:
        n = np.float(len(n))
    except TypeError:
        n = np.float(n)
    if a is None:
        a = 3.0/8 if(n <= 10) else 1.0/2
    return (np.arange(n) + 1 - a)/(n + 1 - 2*a)


def qqnorm(y):
    ina = np.isnan(y)
    if ina.sum() > 0:
        yN = y
        y = y[~ina]
    n = y.shape[0]
    if n == 0:
        print('y is empty or has only NAs')
        return np.array([])
    x = np.around(norm.ppf(ppoints(n)[np.argsort(np.argsort(y))]), decimals=15)
    if ina.sum() > 0:
        y = x
        x = yN
        x[~ina] = y
    return x


def trait_correct(pc, y):
    pc1 = pd.concat([pd.DataFrame(np.ones((y.shape[0], 1)), index=pc.index), pc], axis=1)
    vhat = np.dot(np.linalg.pinv(np.dot(pc1.T, pc1)), np.dot(pc1.T, y))
    if len(vhat.shape) == 1:
        y_corr = y - np.dot(pc, vhat[1:])
    else:
        y_corr = y - np.dot(pc, vhat[1:, :])
    return y_corr


def pc_calc(bed, pc_num):
    try:
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri
        from rpy2.rinterface_lib.embedded import RRuntimeError
        pandas2ri.activate()

        warnings.filterwarnings("ignore")
        base = importr('base')
        bigsnpr = importr('bigsnpr')
        bigstatsr = importr('bigstatsr')

        base.sink('/dev/null')
        g = bigsnpr.snp_readBed(bed, backingfile=base.tempfile())
        g = bigsnpr.snp_attach(g)
        svd = bigsnpr.snp_autoSVD(g[0], infos_chr=g[2]['chromosome'], ncores=1,
                                  infos_pos=g[2]['physical.pos'], thr_r2=np.nan, k=pc_num)
        base.sink()
    except RRuntimeError:
        return None
    else:
        pc = bigstatsr.predict_big_SVD(svd)
        pc = base.data_frame(pc, row_names=g[1]['sample.ID'])
        #pc = base.cbind(pc.rownames, pc)
        pc.columns = ['PC' + str(i) for i in range(1, pc_num+1)]
        return pc