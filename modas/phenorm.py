import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import norm
from sklearn.preprocessing import MinMaxScaler
import warnings
import subprocess
import re


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
        n = np.float64(len(n))
    except TypeError:
        n = np.float64(n)
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
        import rpy2.robjects as robjects
        pandas2ri.activate()

        warnings.filterwarnings("ignore")
        base = importr('base')
        utils = importr('utils')
        if not base.require('bigsnpr', quietly=True)[0]:
            utils_path = subprocess.check_output('locate modas/utils', shell=True, text=True, encoding='utf-8')
            # utils_path = '/'.join(re.search('\n(.*site-packages.*)\n', utils_path).group(1).split('/')[:-1])
            utils_path = re.search('\n(.*site-packages.*)\n', utils_path).group(1)
            if not utils_path.endswith('utils'):
                utils_path = '/'.join(utils_path.split('/')[:-1])
            utils.install_packages(utils_path + '/Matrix_1.6-5.tar.gz', repos=robjects.rinterface.NULL, type='source',
                                   quiet=True)
            utils.install_packages('bigsnpr', dependence=True, repos='https://cloud.r-project.org', quiet=True)
        robjects.r['options'](warn=-1)
        bigsnpr = importr('bigsnpr')
        bigstatsr = importr('bigstatsr')

        base.sink('/dev/null')
        g = bigsnpr.snp_readBed(bed, backingfile=base.tempfile())
        g = bigsnpr.snp_attach(g)
        # svd = bigsnpr.snp_autoSVD(g.rx2[0], infos_chr=g[2]['chromosome'], ncores=1,
        #                           infos_pos=g[2]['physical.pos'], thr_r2=np.nan, k=pc_num)
        svd = bigsnpr.snp_autoSVD(g.rx2('genotypes'), infos_chr=g.rx2('map').rx2('chromosome'),
                                  infos_pos=g.rx2('map').rx2('physical.pos'), thr_r2=np.nan, k=pc_num)
        base.sink()
    except RRuntimeError:
        return None
    else:
        pc = bigstatsr.predict_big_SVD(svd)
        # pc = base.data_frame(pc, row_names=g[1]['sample.ID'])
        #pc = base.cbind(pc.rownames, pc)
        pc = pd.DataFrame(pc, index=g.rx2('fam').rx2('sample.ID'))
        pc.columns = ['PC' + str(i) for i in range(1, pc_num+1)]
        return pc
