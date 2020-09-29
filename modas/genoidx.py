from pandas_plink import read_plink1_bin
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from multiprocessing import cpu_count
import modas.multiprocess as mp
import numpy as np
import pandas as pd
from collections import Counter
import struct




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

def optimal_cluster(g):
    cluster = DBSCAN(eps=0.2, min_samples=10, metric='jaccard').fit(g.T)
    if sum(cluster.labels_ == -1) < g.shape[1]*0.3:
        return cluster
    else:
        cluster = DBSCAN(eps=0.35, min_samples=5, metric='jaccard').fit(g.T)
        return cluster


def cluster_PCA(g,index,colname_prefix):
    g = np.array(g)
    g[g == 2] = 4
    g[g == 0] = 2
    g[g == 4] = 0
    g_std = g - np.mean(g, axis=0)
    cluster = optimal_cluster(g)
    pca = PCA(n_components=3)
    window_res = pd.DataFrame()
    for label in np.unique(cluster.labels_):
        if label != -1 and sum(cluster.labels_ == label) >= 10:
            g_sub = g_std[:, cluster.labels_ == label]
            pca.fit(g_sub)
            pc = pd.DataFrame(pca.transform(g_sub), index=index, columns=[colname_prefix+'_cluster'+str(int(label)+1)+'_PC'+str(i) for i in range(1,4)])
            if pca.explained_variance_ratio_[0] >= 0.6:
                pc = pc.iloc[:,0].to_frame()
            else:
                pc = pc.iloc[:,:2]
            window_res = pd.concat([window_res, pc], axis=1)
    return window_res

def chr_cluster_pca(G_chr,chrom,window,step):
    chr_res = pd.DataFrame()
    paras = []
    chr_end = G_chr.pos[-1]
    start = 0
    end = start+window
    if step == 0:
        step = window
    while start < chr_end:
        G_chr_sub = G_chr.where((G_chr.pos>=start) & (G_chr.pos<=end),drop=True)
        if G_chr_sub.shape[1] <= 100:
            start = start + step
            end = start + window
            continue
        if end > chr_end:
            colname_prefix = '_'.join(['chr',chrom,str(start),str(int(chr_end))])
        else:
            colname_prefix = '_'.join(['chr',chrom,str(start),str(end)])

        paras.append((G_chr_sub, G_chr_sub.sample, colname_prefix))
        cluster_pc = cluster_PCA(G_chr_sub,G_chr_sub.sample,colname_prefix)
        if chr_res.empty:
            chr_res = cluster_pc
        else:
            chr_res = pd.concat([chr_res, cluster_pc], axis=1)
        start = start + step
        end = start + window
    return chr_res


def genome_cluster(G, window, step, threads):
    paras = list()
    if threads > np.unique(G.chrom).shape[0]:
        threads = np.unique(G.chrom).shape[0]
    for chrom in np.unique(G.chrom):
        G_chr = G.where(G.chrom == chrom,drop=True)
        paras.append((G_chr, chrom, window, step))
    res = mp.parallel(chr_cluster_pca, paras, threads)
    res = pd.concat(res, axis=1)
    return res


def read_genotype(geno_prefix):
    try:
        G = read_plink1_bin(geno_prefix + '.bed', geno_prefix + '.bim', geno_prefix + '.fam', verbose=False)
    except Exception:
        return None
    return G


def snp_clumping(bed, r2, out):
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter
    base = importr('base')
    bigsnpr = importr('bigsnpr')
    g = bigsnpr.snp_readBed(bed, backingfile=base.tempfile())
    g = bigsnpr.snp_attach(g)
    with localconverter(ro.default_converter + pandas2ri.converter):
        plink_map = ro.conversion.py2rpy(g[2])
    snp_keep = bigsnpr.snp_clumping(g[0], infos_chr=plink_map.rx2('chromosome'), infos_pos=plink_map.rx2('physical.pos'),
                                    thr_r2=r2)
    g_clump = bigsnpr.subset_bigSNP(g, ind_col=snp_keep)
    g_clump = bigsnpr.snp_attach(g_clump)
    bigsnpr.snp_writeBed(g_clump, out)
    return g[2].nrow, len(snp_keep)
