from pandas_plink import read_plink1_bin
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from multiprocessing import cpu_count
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import pairwise_distances
import warnings
import modas.multiprocess as mp
import numpy as np
import pandas as pd
from collections import Counter
import struct

warnings.filterwarnings("ignore")


def convert(g, n):
    c = [1, 4, 16, 64]
    s = 0
    for i in range(len(g)):
        if g[i][0] != g[i][1]:
            s += 2 * c[i]
        elif g[i] == 'NN':
            s += 1 * c[i]
        elif g[i][0] == n[0]:
            s += 0 * c[i]
        else:
            s += 3 * c[i]
    return s


def hap2plink_bed(hap, plink):
    try:
        with open(hap) as h, open(plink + '.bed', 'wb') as b, open(plink + '.bim', 'w') as bim, open(plink + '.fam','w') as fam:
            #b = open(plink + '.bed', 'wb')
            #bim = open(plink + '.bim', 'w')
            #fam = open(plink + '.fam','w')

            b.write(struct.pack('B', 108))
            b.write(struct.pack('B', 27))
            b.write(struct.pack('B', 1))
            out = list()
            for _, l in enumerate(h):
                l = l.strip().split('\t')
                if _ == 0:
                    samples = l[11:]
                    for s in samples:
                        fam.write(' '.join([s, s, '0', '0', '0', '-9']) + '\n')
                else:
                    g_seq = ''.join(l[11:])
                    nlu_num = list(set(g_seq))
                    if len(nlu_num) > 2 and 'N' not in nlu_num:
                        print("Warning: There is more than two nucleotide letter snp in hapmap file. skip this snp")
                        continue
                    if len(nlu_num) == 1:
                        print("Warning: There is one nucleotide letter snp in hapmap file. skip this snp")
                        continue
                    if 'N' in nlu_num:
                        c = Counter(g_seq)
                        del c['N']
                        n = sorted(c, key=lambda x: c[x])
                    else:
                        g_seq_sort = ''.join(sorted(g_seq))
                        major = g_seq_sort[len(g_seq) // 2]
                        if g_seq_sort[0] == major:
                            n = [g_seq_sort[-1], major]
                        else:
                            n = [g_seq_sort[0], major]
                    bim.write('\t'.join([l[2], l[0], '0', l[3]] + n) + '\n')
                    for num in range(11, len(l), 4):
                        if num + 4 < len(l):
                            #out += struct.pack('B', convert(l[num:num + 4], n)).decode()
                            out.append(convert(l[num:num + 4], n))
                        else:
                            #out += struct.pack('B', convert(l[num:len(l)], n)).decode()
                            out.append(convert(l[num:len(l)], n))
            #b.write(out)
            b.write(struct.pack('B'*len(out),*out))
    except Exception:
        return False
    else:
        return True


#def jaccard_tri_new(u, v):
#    a = (u == v).sum() - np.bitwise_and(u == 0, v == 0).sum()
#    b = np.logical_or(u, v).sum()
#    return 1 - a / b if b != 0 else 0


def jaccard_tri(u, v):
    u = u.astype(int)
    v = v.astype(int)
    a = np.logical_and(u, v)
    minor_het_u = (np.bitwise_and(u, v)).astype(bool).sum()
    minor_het_i = (np.bitwise_or(u, v)).astype(bool).sum() + a.sum() - minor_het_u
    return 1 - np.double(minor_het_u)/minor_het_i if minor_het_i !=0 else 0



def optimal_cluster(g):
    dis = pairwise_distances(g.T, metric=jaccard_tri, n_jobs=-1)
    cluster = DBSCAN(eps=0.2, min_samples=10, metric='precomputed').fit(dis)
    if sum(cluster.labels_ == -1) < g.shape[1]*0.3:
        return cluster
    else:
        # cluster = DBSCAN(eps=0.35, min_samples=5, metric='jaccard').fit(dis)
        #cluster = DBSCAN(eps=0.35, min_samples=5, metric=jaccard_tri_new).fit(g.T)
        # cluster = DBSCAN(eps=0.35, min_samples=5, metric='jaccard').fit(g.T)
        cluster = DBSCAN(eps=0.35, min_samples=5, metric='precomputed').fit(dis)
        return cluster


def cluster_PCA(g,index,colname_prefix):
    g = g.values
    g_std = g - np.mean(g, axis=0)
    cluster = optimal_cluster(g)
    pca = PCA(n_components=3)
    window_res = pd.DataFrame()
    #cluster_variant = pd.DataFrame()
    for label in np.unique(cluster.labels_):
        if label != -1 and sum(cluster.labels_ == label) >= 10:
            g_sub = g_std[:, cluster.labels_ == label]
            try:
                pca.fit(g_sub)
            except np.linalg.LinAlgError:
                continue
            pc = pd.DataFrame(pca.transform(g_sub), index=index, columns=[colname_prefix+'_cluster'+str(int(label)+1)+'_PC'+str(i) for i in range(1,4)])
            if pca.explained_variance_ratio_[0] >= 0.6:
                pc = pc.iloc[:,0].to_frame()
            else:
                pc = pc.iloc[:,:2]
            window_res = pd.concat([window_res, pc], axis=1)
            #variant = pd.Series(g.snp.values)[cluster.labels_ == label].to_frame()
            #variant.loc[:,'cluster'] = colname_prefix+'_cluster'+str(int(label)+1)
            #variant.columns = ['rs','cluster']
            #cluster_variant = pd.concat([cluster_variant,variant])
    #return window_res,cluster_variant
    return window_res

def chr_cluster_pca(G_chr,chrom,window,step):
    chr_res = pd.DataFrame()
    #chr_cluster_variant = pd.DataFrame()
    paras = []
    chr_end = G_chr.pos[-1]
    start = 0
    end = start+window
    if step == 0:
        step = window
    while start < chr_end:
        G_chr_sub = G_chr.where((G_chr.pos>=start) & (G_chr.pos<=end), drop=True)
        if G_chr_sub.shape[1] <= 100:
            start = start + step
            end = start + window
            continue
        if end > chr_end:
            colname_prefix = '_'.join(['chr',chrom,str(start),str(int(chr_end))])
        else:
            colname_prefix = '_'.join(['chr',chrom,str(start),str(end)])

        paras.append((G_chr_sub, G_chr_sub.sample, colname_prefix))
        #cluster_pc,cluster_variant = cluster_PCA(G_chr_sub,G_chr_sub.sample,colname_prefix)
        cluster_pc = cluster_PCA(G_chr_sub,G_chr_sub.sample,colname_prefix)
        if chr_res.empty:
            chr_res = cluster_pc
        else:
            chr_res = pd.concat([chr_res, cluster_pc], axis=1)
        #chr_cluster_variant = pd.concat([chr_cluster_variant,cluster_variant])
        start = start + step
        end = start + window
    #return chr_res,chr_cluster_variant
    return chr_res


def genome_cluster(G, window, step, threads):
    paras = list()
    if threads > np.unique(G.chrom).shape[0]:
        threads = np.unique(G.chrom).shape[0]
    for chrom in np.unique(G.chrom):
        G_chr = G.where(G.chrom == chrom,drop=True)
        paras.append((G_chr, chrom, window, step))
    res = mp.parallel(chr_cluster_pca, paras, threads)
    #res_pc = pd.concat([i[0] for i in res], axis=1)
    res_pc = pd.concat(res, axis=1)
    res_pc.loc[:,:] = np.around(MinMaxScaler(feature_range=(0, 2)).fit_transform(res_pc.values),decimals=3)
    #res_variant = pd.concat([i[1] for i in res])
    #return res_pc,res_variant
    return res_pc


def read_genotype(geno_prefix):
    try:
        G = read_plink1_bin(geno_prefix + '.bed', geno_prefix + '.bim', geno_prefix + '.fam', ref='a0',verbose=False)
    except Exception:
        return None
    return G

def snp_clumping(bed, r2, out):
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    import rpy2.robjects as robjects

    pandas2ri.activate()
    robjects.r['options'](warn=-1)
    robjects.r('options(datatable.showProgress = FALSE)')
    base = importr('base')
    base.sink('/dev/null')
    bigsnpr = importr('bigsnpr')
    g = bigsnpr.snp_readBed(bed, backingfile=base.tempfile())
    g = bigsnpr.snp_attach(g)
    snp_keep = bigsnpr.snp_clumping(g[0], infos_chr=g[2]['chromosome'], infos_pos=g[2]['physical.pos'],
                                    thr_r2=r2, ncores=1)
    g_clump = bigsnpr.subset_bigSNP(g, ind_col=snp_keep)
    g_clump = bigsnpr.snp_attach(g_clump)
    bigsnpr.snp_writeBed(g_clump, out)
    base.sink()
    return g[2].shape[0], len(snp_keep)
