import pandas as pd
import numpy as np
import bioframe as bf
from image_match.goldberg import ImageSignature
from pandas_plink import read_plink1_bin
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, leaves_list, cut_tree
from joblib import Parallel, delayed
from sklearn.metrics import silhouette_score, calinski_harabasz_score
import resource
import os

resource.setrlimit(resource.RLIMIT_NOFILE, (4096, 4096))


def qtl_cluster(qtl):
    qtl.CHR = qtl.CHR.astype(str)
    qtl = bf.cluster(qtl, cols=['CHR', 'qtl_start', 'qtl_end'], min_dist=0)
    return qtl


def kin(g, top=2):
    g = g - g.mean()
    K = np.dot(g, g.T)
    d = np.diag(K)
    DL = np.min(d)
    DU = np.max(d)
    floor = np.min(K)
    K = top * (K - floor) / (DU - floor)
    Dmin = top * (DL - floor) / (DU - floor)
    dig_index = np.eye(K.shape[0], dtype=bool)
    if Dmin < 1:
        K[dig_index] = (np.diag(K)-Dmin+1)/((top+1-Dmin)*0.5)
        K[~dig_index] = K[~dig_index] * (1 / Dmin)
    Omax = np.max(K[~dig_index])
    if Omax > top:
        K[~dig_index] = K[~dig_index] * (top / Omax)
    return K


def get_kin_info(qtl, gwas_dir, geno, pvalue):
    var = pd.Series(geno.variant, index=geno.snp)
    kin_info = dict()
    for phe_name in qtl.phe_name.unique():
        fn = gwas_dir+'/tmp_' + phe_name + '_plink.assoc.txt'
        if not os.path.exists(fn):
            print('Warning: ' + fn + 'is not exist.')
            continue
        gwas = pd.read_csv(fn, sep='\t')
        geno_sub = geno.sel(variant=var.reindex(gwas.loc[gwas.p_wald <= pvalue, 'rs']).dropna().values, drop=True)
        geno_sub = pd.DataFrame(geno_sub.values, index=geno_sub.fid, columns=geno_sub.snp)
        kin_res = kin(geno_sub)
        ril_cluster = linkage(kin_res, method='ward')
        idx = leaves_list(ril_cluster)
        label = cut_tree(ril_cluster, n_clusters=2)[:, 0]
        kin_info[phe_name] = dict([['kin', kin_res], ['idx', idx], ['label', label]])
    return kin_info


def get_ril_cluster_idx(kin1, kin2, metric):
    score1 = calc_cluster_score(kin1['kin'], kin2['kin'], kin1['label'], metric)
    score2 = calc_cluster_score(kin1['kin'], kin2['kin'], kin2['label'], metric)
    if metric == 'silhouette':
        if score1 > score2:
            return kin1['idx']
        else:
            return kin2['idx']
    if metric == 'calinski_harabasz':
        if score1 > score2:
            return kin1['idx']
        else:
            return kin2['idx']


def calc_cluster_score(kin1, kin2, label, metric):
    if metric == 'silhouette':
        from sklearn.metrics import silhouette_score
        kin1_score = silhouette_score(kin1, label)
        kin2_score = silhouette_score(kin2, label)
    elif metric == 'calinski_harabasz':
        from sklearn.metrics import calinski_harabasz_score
        kin1_score = calinski_harabasz_score(kin1, label)
        kin2_score = calinski_harabasz_score(kin2, label)
    return np.mean([kin1_score, kin2_score])


def get_signature(g, gis):
    # image_limits = gis.crop_image(g,
    #                               lower_percentile=gis.lower_percentile,
    #                               upper_percentile=gis.upper_percentile,
    #                               fix_ratio=gis.fix_ratio)
    # x_coords, y_coords = gis.compute_grid_points(g,
    #                                              n=gis.n, window=image_limits)
    x_coords, y_coords = gis.compute_grid_points(g, n=gis.n)
    avg_grey = gis.compute_mean_level(g, x_coords, y_coords, P=gis.P)
    diff_mat = gis.compute_differentials(avg_grey,
                                         diagonal_neighbors=gis.diagonal_neighbors)
    gis.normalize_and_threshold(diff_mat,
                                identical_tolerance=gis.identical_tolerance,
                                n_levels=gis.n_levels)
    return np.ravel(diff_mat).astype('int8')


def calc_image_match_score(kin_info, phe_list, metric):
    gis = ImageSignature()
    score = list()
    for _, phe1 in enumerate(phe_list):
        for phe2 in phe_list[_+1:]:
            idx = get_ril_cluster_idx(kin_info[phe1], kin_info[phe2], metric)
            score.append(gis.normalized_distance(get_signature(kin_info[phe1]['kin'][idx, :][:, idx], gis), get_signature(kin_info[phe2]['kin'][idx, :][:, idx], gis)))
    score = squareform(score)
    return score


def cluster_coloc(kin_info, qtl, c, metric, cls):
    qtl_sub = qtl.loc[qtl.cluster==c, :]
    if qtl_sub.shape[0] > 1:
        phe_name = qtl_sub.phe_name.unique()
        trait_dis = calc_image_match_score(kin_info, phe_name, metric)
        trait_dis = np.round(trait_dis, 2)
        cls.fit(trait_dis)
        cls_res = pd.Series(cls.labels_, index=phe_name).to_frame().reset_index()
        cls_res.columns = ['phe_name', 'label']
        cls_res = cls_res.loc[cls_res.label != -1, :]
        if not cls_res.empty:
            cls_res['label'] = str(c) + '_' + cls_res.label.astype(str)
            return pd.DataFrame(trait_dis, index=phe_name, columns=phe_name), cls_res
        else:
            return pd.DataFrame(trait_dis, index=phe_name, columns=phe_name), pd.DataFrame()
    else:
        return pd.DataFrame(), pd.DataFrame()


def trait_coloc(kin_info, qtl, metric, eps, p):
    cls = DBSCAN(eps=eps, min_samples=2, metric='precomputed')
    # dis = list()
    # coloc_res = list()
    # coloc_count = 0
    # for c in qtl.cluster.unique():
    #     qtl_sub = qtl.loc[qtl.cluster==c, :]
    #     if qtl_sub.shape[0] > 1:
    #         phe_name = qtl_sub.phe_name.unique()
    #         trait_dis = calc_image_match_score(kin_info, phe_name, metric)
    #         dis.append(pd.DataFrame(trait_dis, index=phe_name, columns=phe_name))
    #         cls.fit(trait_dis)
    #         cls_res = pd.Series(cls.labels_, index=phe_name).to_frame().reset_index()
    #         cls_res.columns = ['phe_name', 'label']
    #         cls_res = cls_res.loc[cls_res.label != -1, :]
    #         if not cls_res.empty:
    #             cls_count = cls_res['label'].value_counts()
    #             cls_res['label'] = cls_res['label'].replace(cls_count.index, np.arange(coloc_count + 1, coloc_count + 1 + cls_count.shape[0]))
    #             coloc_res.append(cls_res)
    #             coloc_count = coloc_count + cls_count.shape[0]
    res = Parallel(n_jobs=p)(delayed(cluster_coloc)(kin_info, qtl, c, metric, cls) for c in qtl.cluster.unique())
    coloc_res = [i[1] for i in res]
    dis = [i[0] for i in res]
    coloc_res = pd.concat(coloc_res, axis=0)
    if not coloc_res.empty:
        coloc_res_count = coloc_res['label'].value_counts()
        coloc_res['label'] = coloc_res['label'].replace(coloc_res_count.index, np.arange(1, coloc_res_count.shape[0] + 1))
        coloc_res = pd.merge(qtl.drop(['cluster', 'cluster_start', 'cluster_end'], axis=1), coloc_res, on='phe_name', how='left')
        coloc_res = coloc_res.fillna(-1)
        coloc_res = coloc_res.loc[coloc_res.label != -1, :]
    dis = pd.concat(dis)
    dup_index = dis.index[dis.index.duplicated()]
    dup_dis = dis[dis.index.duplicated(keep=False)]
    dis = dis[~dis.index.duplicated()]
    for index in dup_index:
        dis.loc[index, :] = dup_dis.loc[index, :].apply(lambda x:pd.Series(x[~pd.isna(x)]).min() if(pd.isna(x).sum()!=x.shape[0]) else x[0], axis=0)
    dis = dis.fillna(1)
    dis_pairwise = dis.stack().reset_index()
    dis_pairwise.columns = ['level_0', 'level_1', 'image_match_score']
    dis_pairwise = dis_pairwise.loc[dis_pairwise.level_0 != dis_pairwise.level_1, :]
    dis_pairwise['id'] = dis_pairwise.apply(lambda x: ';'.join(sorted([x['level_0'], x['level_1']])), axis=1)
    dis_pairwise = dis_pairwise.drop_duplicates(subset='id')
    qtl_overlap = bf.overlap(qtl[['CHR', 'qtl_start', 'qtl_end', 'SNP', 'P', 'phe_name']], qtl[['CHR', 'qtl_start', 'qtl_end', 'SNP', 'P', 'phe_name']],
                             cols1=['CHR', 'qtl_start', 'qtl_end'], cols2=['CHR', 'qtl_start', 'qtl_end'], how='inner', suffixes=('_1', '_2'))
    qtl_overlap = qtl_overlap.loc[qtl_overlap.phe_name_1 != qtl_overlap.phe_name_2, :]
    qtl_overlap['id'] = qtl_overlap.apply(lambda x: ';'.join(sorted([x['phe_name_1'], x['phe_name_2']])), axis=1)
    qtl_overlap = qtl_overlap.drop_duplicates(subset='id')
    coloc_pairwise_res = pd.merge(qtl_overlap, dis_pairwise, on='id')
    coloc_pairwise_res = coloc_pairwise_res.drop(['id', 'level_0', 'level_1'], axis=1)
    coloc_pairwise_res['coloc'] = 'No'
    coloc_pairwise_res.loc[coloc_pairwise_res['image_match_score'] <= 0.2, 'coloc'] = 'Yes'
    dis = dis.reindex(qtl.phe_name.unique()).reindex(qtl.phe_name.unique(), axis=1)
    dis = dis.fillna(1)
    return coloc_res, coloc_pairwise_res, dis

