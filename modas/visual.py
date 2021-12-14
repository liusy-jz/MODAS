import pandas as pd
import numpy as np
import modas.multiprocess as mp
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.embedded import RRuntimeError
import rpy2.robjects as robjects
from collections import Counter
import modas.gwas_cmd as gc
import pyranges as pr
from yattag import Doc, indent
import shutil
import warnings
import glob
import os
import re

pandas2ri.activate()
data_table = importr('data.table')
base = importr('base')
robjects.r('options(datatable.showProgress = FALSE)')

warnings.filterwarnings("ignore")

# def gwas(phe, geno, num_threads, phe_fn):
#     geno_prefix = geno.split('/')[-1]
#     related_matrix_cmd = 'gemma.linux -bfile {0}.link -gk 1 -o {1}'.format(geno_prefix,geno_prefix)
#     gwas_cmd = 'gemma.linux -bfile {0}.link -k output/{0}.cXX.txt -lmm -n {1} -o {2}'
#     fam = pd.read_csv(geno+'.fam', sep=r'\s+', header=None)
#     fam[5] = 1
#     fam = pd.merge(fam, phe, left_on=0, right_index=True, how='left')
#     fam.to_csv(geno_prefix+'.link.fam', sep='\t', na_rep='NA', header=None, index=False)
#     if os.path.exists(geno_prefix+'.link.bed'):
#         os.remove(geno_prefix+'.link.bed')
#     if os.path.exists(geno_prefix+'.link.bim'):
#         os.remove(geno_prefix+'.link.bim')
#     os.symlink(geno+'.bed', geno_prefix+'.link.bed')
#     os.symlink(geno+'.bim', geno_prefix+'.link.bim')
#     values = list()
#     for _, p in enumerate(phe.columns):
#         p = p.replace('/', '.')
#         values.append((gwas_cmd.format(*[geno_prefix, _ + 2, '.'.join(phe_fn.split('/')[-1].split('.')[:-1])+ '_' + str(p)]),))
#     s = mp.run(related_matrix_cmd)
#     if s != 0:
#         return None
#     else:
#         s = mp.parallel(mp.run, values, num_threads)
#         os.remove(geno_prefix+'.link.bed')
#         os.remove(geno_prefix+'.link.bim')
#         os.remove(geno_prefix+'.link.fam')
#         return s


def gwas(phe, geno, num_threads, phe_fn, gwas_model):
    software, model = gwas_model.split('_')
    geno_prefix = geno.split('/')[-1]
    phe_fn = '.'.join(phe_fn.split('/')[-1].split('.')[:-1])
    if software == 'gemma' and model == 'MLM':
        geno_prefix = geno.split('/')[-1]
        related_matrix_cmd = 'gemma.linux -bfile {0}.link -gk 1 -o {1}'.format(geno_prefix,geno_prefix)
        fam = pd.read_csv(geno + '.fam', sep=r'\s+', header=None)
        fam[5] = 1
        fam = pd.merge(fam, phe, left_on=0, right_index=True, how='left')
        fam.to_csv(geno_prefix+'.link.fam', sep='\t', na_rep='NA', header=None, index=False)
    if software != 'GAPIT' and gwas_model != 'gemma_LM':
        if os.path.exists(geno_prefix+'.link.bed'):
            os.remove(geno_prefix+'.link.bed')
        if os.path.exists(geno_prefix+'.link.bim'):
            os.remove(geno_prefix+'.link.bim')
        os.symlink(geno+'.bed', geno_prefix+'.link.bed')
        os.symlink(geno+'.bim', geno_prefix+'.link.bim')
    if software == 'rMVP':
        if os.path.exists(geno_prefix+'.link.fam'):
            os.remove(geno_prefix+'.link.fam')
        os.symlink(geno + '.fam', geno_prefix + '.link.fam')
        fam = pd.read_csv(geno + '.fam', sep=r'\s+', header=None)
        phe = phe.reindex(fam[0])
    if software == 'gemma' and model == "LM":
        if os.path.exists('./gemma_lm_geno'):
            shutil.rmtree('./gemma_lm_geno')
        os.mkdir('./gemma_lm_geno')
        for p in phe.columns:
            p = p.replace('m/z', 'm.z')
            os.symlink('../' + geno + '.bed', './gemma_lm_geno/' + geno_prefix + '_' + p + '.link.bed')
            os.symlink('../' + geno + '.bim', './gemma_lm_geno/' + geno_prefix + '_' + p + '.link.bim')
            fam = pd.read_csv(geno + '.fam', sep=r'\s+', header=None)
            fam = pd.merge(fam.iloc[:, :5], phe.loc[:, p.replace('m.z', 'm/z')].to_frame(), left_on=0, right_index=True, how='left')
            fam.to_csv('./gemma_lm_geno/' + geno_prefix + '_' + p + '.link.fam', sep='\t', na_rep='NA', header=None, index=False)

    values = list()
    for _, p in enumerate(phe.columns):
        p = p.replace('m/z', 'm.z')
        if gwas_model == 'gemma_LM':
            values.append(((gc.gemma_cmd(model, './gemma_lm_geno/' + geno_prefix + '_' + p + '.link', None, None, '_'.join([phe_fn, gwas_model, p]))),))
        if gwas_model == 'gemma_MLM':
            values.append(((gc.gemma_cmd(model, geno_prefix + '.link', geno_prefix, _ + 2, '_'.join([phe_fn, gwas_model, p]))),))
        if software == 'rMVP':
            if not os.path.exists('./output'):
                os.mkdir('./output')
            omics_phe = phe.loc[:, p.replace('m.z', 'm/z')].to_frame().reset_index()
            omics_phe.columns = ['Taxa', '_'.join([phe_fn, gwas_model, p])]
            values.append((model, geno_prefix+'.link', geno_prefix+'.link', omics_phe, 1, './output'))
    if gwas_model == 'gemma_MLM':
        s = mp.run(related_matrix_cmd)
        if s != 0:
            return None
    if software == 'gemma':
        s = mp.parallel(mp.run, values, num_threads)
    if software == 'rMVP':
        s1 = mp.parallel(gc.rmvp, (values[0],), 1)
        s = mp.parallel(gc.rmvp, values[1:], num_threads)
        s = s1 + s
    if software == 'GAPIT':
        if not os.path.exists('./output'):
            os.mkdir('./output')
        for path in os.environ.get('PATH').split(':'):
            if re.search(r'MODAS/utils', path):
                gapit_path = path
        phe.columns = ['_'.join([phe_fn, gwas_model, p]) for p in phe.columns]
        phe = phe.reset_index()
        phe.columns = ['Taxa'] + list(phe.columns[1:])
        geno = os.path.abspath(geno)
        os.chdir('./output')
        s = gc.gapit(model, geno, phe, gapit_path)
        os.chdir('../')
    if gwas_model == 'gemma_MLM' or software == 'rMVP':
        os.remove(geno_prefix+'.link.bed')
        os.remove(geno_prefix+'.link.bim')
        os.remove(geno_prefix+'.link.fam')
    if gwas_model == 'gemma_LM':
        shutil.rmtree('./gemma_lm_geno')
    return s


def gwas_plot(res, p, prefix, t, software):
    try:
        base.sink('/dev/null')
        w = data_table.fread(res, data_table=base.getOption("datatable.fread.datatable", False))
        if software == 'gemma':
            w_subset = w.loc[w.p_wald <= float(p), :]
            m = w_subset[['rs', 'chr', 'ps', 'p_wald']]
            q = w[['rs', 'chr', 'ps', 'p_wald']]
        if software == 'rMVP':
            w_subset = w.loc[w[w.columns[-1]] <= float(p), :]
            m = w_subset.iloc[:, [0, 1, 2, -1]]
            q = w.iloc[:, [0, 1, 2, -1]]
        if software == 'GAPIT':
            w_subset = w.loc[w['P.value'] <= float(p), :]
            m = w_subset[['SNP', 'Chromosome', 'Position', 'P.value']]
            q = w[['SNP', 'Chromosome', 'Position', 'P.value']]
        m.columns = ['SNP', 'Chromosome', 'Position', prefix]
        q.columns = ['SNP', 'Chromosome', 'Position', prefix]
        thresholdi = robjects.FloatVector([1.0 / w.shape[0], 1e-6, 1e-5])
        lim = -np.log10(min(m[prefix])) + 2
        #w_subset = base.subset(w, np.array(w.rx2('p_wald')) <= float(p))
        #m = w_subset.rx(robjects.StrVector(['rs', 'chr', 'ps', 'p_wald']))
        #q = w.rx(robjects.StrVector(['rs', 'chr', 'ps', 'p_wald']))
        # m.names = ['SNP', 'Chromosome', 'Position', prefix]
        # q.names = ['SNP', 'Chromosome', 'Position', prefix]
        #thresholdi = robjects.FloatVector([1.0/w.nrow, 1e-6, 1e-5])
        #lim = -np.log10(min(np.array(w_subset.rx2('p_wald'))))+2
        #base.sink('/dev/null')
        for path in os.environ.get('PATH').split(':'):
            if re.search(r'MODAS/utils',path):
                robjects.r('source("'+path+'/CMplot.r")')
        CMplot = robjects.r['CMplot']
        CMplot(m, plot_type='m', col=robjects.StrVector(["grey30", "grey60"]), ylim=robjects.FloatVector([2, lim]), threshold=thresholdi,
                cex=robjects.FloatVector([0.5, 0.5, 0.5]), signal_cex=robjects.FloatVector([0.5, 0.5, 0.5]),
                threshold_col=robjects.StrVector(['red', 'green', 'blue']), chr_den_col=robjects.rinterface.NULL, amplify=True,
                signal_pch = robjects.IntVector([19, 19, 19]), dpi=300,
                signal_col=robjects.StrVector(['red', 'green', 'blue']), multracks=False, LOG10=True, file=t)
        CMplot(q, plot_type='q', col='grey30', threshold=thresholdi[0],
               signal_cex=robjects.FloatVector([0.5, 0.5, 0.5]), signal_pch=robjects.IntVector([19, 19, 19]),
               conf_int_col='gray', signal_col='red', multracks=False, LOG10=True, file=t, dpi=300)
        base.sink()
    except RRuntimeError:
        return 0
    except ValueError:
        return 0
    else:
        return 1


def gwas_plot_parallel(phe, p, threads, t, phe_fn, gwas_model):
    software, model = gwas_model.split('_')
    values = list()
    for i in phe.columns:
        i = str(i).replace('m/z', 'm.z')
        if software == 'gemma':
            gwas_fn = 'output/' + '_'.join(['.'.join(phe_fn.split('/')[-1].split('.')[:-1]), gwas_model, str(i)]) + '.assoc.txt'
        if software == 'rMVP':
            gwas_fn = 'output/' + '_'.join(['.'.join(phe_fn.split('/')[-1].split('.')[:-1]), gwas_model, str(i)]) + '.'.join(['.'+model, 'csv'])
        if software == 'GAPIT':
            gwas_fn = 'output/' + '.'.join(['GAPIT', model, '_'.join(['.'.join(phe_fn.split('/')[-1].split('.')[:-1]), gwas_model, str(i)]), 'GWAS', 'Results', 'csv'])
        values.append((gwas_fn, p, '.'.join(phe_fn.split('/')[-1].split('.')[:-1]) + '_' + str(i), t, software))
    s = mp.parallel(gwas_plot, values, threads)
    return s


def boxplot(phe, g, qtl):
    robjects.r('''box_plot <- function(d, phe, rs, level){
    library(ggplot2)
    library(ggsignif)
    d <- d[d$haplotype!=1,]
    d[d$haplotype==0,'haplotype'] <- level[1]
    d[d$haplotype==2,'haplotype'] <- level[2]
    d$haplotype <- factor(d$haplotype,levels = level)
    b <- as.numeric(formatC(max(d[,1],na.rm=T)*1.2/4,format = 'e',digits = 1))
    p <- ggplot(data = d,aes_string(x='haplotype',y=names(d)[1],fill='haplotype'))+
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size = 4),
          legend.key.size = unit(3,'mm'),
          legend.position = 'none',
          plot.title = element_text(hjust=0.5,size = 6),
          plot.margin=unit(c(0.3,0.3,0,0),'cm'),
          panel.grid = element_blank(),
          axis.line = element_line(colour = 'black',size=0.4),
          axis.text = element_text(size = 6,color = 'black'),
          axis.ticks.length=unit(.1, 'cm'))+
    stat_boxplot(geom = 'errorbar', width = 0.2,size=0.1)+
    geom_boxplot(lwd=0.2,width=0.5,outlier.size = 0.2)+
    geom_signif(comparisons = list(level),map_signif_level = F,
                test= t.test, size=0.2 ,textsize=2, y_position = max(d[,1],na.rm=T)*1.1)+
    #xlab('')+ylab('')+ggtitle(paste(phe,rs,sep='_'))+
    xlab('')+ylab('')+ggtitle('')+
    scale_y_continuous(breaks=seq(0,4*b,by=b),labels = function(x) formatC(x, format = 'e',digits = 1), limits = c(0, max(d[,1], na.rm=T)*1.2))+
    scale_fill_manual(values=c('#E3FFE2', 'forest green'))
    ggsave(paste(phe,'_',rs,'_','boxplot','.jpg',sep=''),plot=p,device='jpg',width=3.5,height=4.3,units = 'cm')
}''')
    robjects.r['options'](warn=-1)
    base.sink('/dev/null')
    box_plot = robjects.r['box_plot']
    g = g.where(g.snp.isin(qtl.SNP), drop=True)
    allele = pd.DataFrame([g.a0, g.a1], index=['a0', 'a1'], columns=g.snp.values)
    g = pd.DataFrame(g.values, index=g.sample, columns=g.snp.values)
    ril = g.index.intersection(phe.index)
    g = g.reindex(ril)
    phe = phe.reindex(ril)
    for index, row in qtl.iterrows():
        if row['phe_name'] not in phe.columns:
            continue
        d = pd.concat([phe[row['phe_name']], g[row['SNP']]], axis=1)
        d.columns = ['trait.' + d.columns[0].replace('-', '.'), 'haplotype']
        level = robjects.StrVector([allele[row['SNP']]['a1'].values*2, allele[row['SNP']]['a0'].values*2])
        box_plot(d, row['phe_name'], row['SNP'], level)
    base.sink()


def multi_trait_plot(phe, gwas_dir, qtl, phe_fn, prefix, t, gwas_model):
    software, model = gwas_model.split('_')
    bk = pd.DataFrame()
    for i in phe.columns:
        i = i.replace('/', '.')
        # fn = gwas_dir + '/' + '.'.join(phe_fn.split('/')[-1].split('.')[:-1]) + '_' + str(i)+'.assoc.txt'
        if software == 'gemma':
            gwas_fn = gwas_dir + '_'.join(['.'.join(phe_fn.split('/')[-1].split('.')[:-1]), gwas_model, str(i)]) + '.assoc.txt'
            d = pd.read_csv(gwas_fn, sep='\t')
            d = d[['rs', 'chr', 'ps', 'p_wald']]
        if software == 'rMVP':
            gwas_fn = gwas_dir + '_'.join(['.'.join(phe_fn.split('/')[-1].split('.')[:-1]), gwas_model, str(i)]) + '.'.join(['.' + model, 'csv'])
            d = pd.read_csv(gwas_fn)
            d = d.iloc[:, [0, 1, 2, -1]]
            d.columns = ['rs', 'chr', 'ps', 'p_wald']
        if software == 'GAPIT':
            gwas_fn = gwas_dir + '.'.join(['GAPIT', model, '_'.join(['.'.join(phe_fn.split('/')[-1].split('.')[:-1]), gwas_model, str(i)]), 'GWAS', 'Results', 'csv'])
            d = pd.read_csv(gwas_fn)
            d = d[['SNP', 'Chromosome', 'Position ', 'P.value']]
            d.columns = ['rs', 'chr', 'ps', 'p_wald']

        # d = pd.read_csv(gwas_fn, sep='\t')
        if bk.empty:
            bk = d.copy()
            bk.loc[bk.p_wald <= 1e-5, 'p_wald'] = 1e-5
        for index, row in qtl.loc[qtl.phe_name == i, :].iterrows():
            peak_pos = int(row['SNP'].split('_')[-1])
            chrom = row['CHR']
            sig_tmp = pd.concat([bk.loc[(bk.chr.astype(str) == str(chrom)) & (bk.ps >= peak_pos-1000000) & (bk.ps <= peak_pos+1000000), 'p_wald'],
                                 d.loc[(bk.chr.astype(str) == str(chrom)) & (bk.ps >= peak_pos-1000000) & (bk.ps <= peak_pos+1000000), 'p_wald']], axis=1)
            sig_tmp.columns = ['bk', 'phe']
            bk.loc[(bk.chr.astype(str) == str(chrom)) & (bk.ps >= peak_pos-1000000) & (bk.ps <= peak_pos+1000000), 'p_wald'] = sig_tmp.apply(
                lambda x: x['bk'] if x['bk'] < x['phe'] else x['phe'], axis=1)
    bk = bk.loc[bk.p_wald <= 1e-2, :]
    bk.loc[bk.p_wald <= 1e-20, 'p_wald'] = 1e-20
    bk = bk[['rs', 'chr', 'ps', 'p_wald']]
    thresholdi = robjects.FloatVector([1.0 / d.shape[0], 1e-6, 1e-5])
    lim = -np.log10(min(bk['p_wald'])) + 2
    bk.columns = ['SNP', 'Chromosome', 'Position', prefix]
    base.sink('/dev/null')
    for path in os.environ.get('PATH').split(':'):
        if re.search(r'MODAS/utils', path):
            robjects.r('source("'+path+'/CMplot.r")')
    #robjects.r('source("/home/debian/文档/MGWAP/compound_extract/plot/CMplot.r")')
    CMplot = robjects.r['CMplot']
    CMplot(bk, plot_type='m', col=robjects.StrVector(["grey30", "grey60"]), ylim=robjects.FloatVector([2, lim]),
           threshold=thresholdi,
           cex=robjects.FloatVector([0.5, 0.5, 0.5]), signal_cex=robjects.FloatVector([0.5, 0.5, 0.5]),
           threshold_col=robjects.StrVector(['red', 'green', 'blue']), chr_den_col=robjects.rinterface.NULL,
           amplify=True,
           signal_pch=robjects.IntVector([19, 19, 19]), dpi=300,
           signal_col=robjects.StrVector(['red', 'green', 'blue']), multracks=False, LOG10=True, file=t)
    base.sink()


def qtl_anno(qtl, anno):
    anno = anno[['geneid', 'position']]
    anno.loc[:, 'chr'] = anno['position'].apply(lambda x: x.split(':')[0])
    anno.loc[:, 'start'] = anno['position'].apply(lambda x: x.split(':')[1].split('-')[0])
    anno.loc[:, 'end'] = anno['position'].apply(lambda x: x.split(':')[1].split('-')[1])
    anno = anno[['chr', 'start', 'end', 'geneid']]
    anno.columns = ['Chromosome', 'Start', 'End', 'geneid']
    qtl_range = qtl[['CHR', 'qtl_start', 'qtl_end', 'phe_name']]
    qtl_range.columns = ['Chromosome', 'Start', 'End', 'phe_name']
    qtl_range = pr.PyRanges(qtl_range)
    anno = pr.PyRanges(anno)
    qtl_anno_intersect = qtl_range.join(anno, how='left')
    qtl_anno = pd.DataFrame()
    for k in sorted(qtl_anno_intersect.dfs.keys()):
        qtl_anno = pd.concat([qtl_anno, qtl_anno_intersect.dfs[k]])
    qtl_anno = qtl_anno[['Chromosome', 'Start', 'End', 'phe_name', 'geneid']]
    qtl_anno.loc[:, 'Chromosome'] = qtl_anno.Chromosome.astype(str)
    qtl_anno.loc[:, 'Start'] = qtl_anno.Start.astype(int)
    qtl_anno.loc[:, 'End'] = qtl_anno.End.astype(int)
    qtl.loc[:, 'CHR'] = qtl.CHR.astype(str)
    qtl = pd.merge(qtl, qtl_anno, left_on=['CHR', 'qtl_start', 'qtl_end', 'phe_name'],
                   right_on=['Chromosome', 'Start', 'End', 'phe_name'])
    qtl = qtl.drop(['Chromosome', 'Start', 'End'], axis=1)
    qtl = qtl.groupby(['CHR', 'qtl_start', 'qtl_end', 'SNP', 'P', 'SP2_num', 'qtl_length', 'phe_name'])['geneid'].apply(';'.join).reset_index()
    qtl.columns = ['CHR', 'qtl_start', 'qtl_end', 'SNP', 'P', 'SP2_num', 'qtl_length', 'phe_name', 'qtl_all_gene']
    qtl.loc[qtl.qtl_all_gene == '-1', 'qtl_all_gene'] = 'nan'
    return qtl


'''
Generate html report for SingleTrait
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2020-08-26 20:33:52
Last modified: 2021-03-24 17:20:33
'''


################# Classes #################
class AllQtlStatistics():
    def __init__(self):
        self.totalSNPs = 0
        self.qtlDetected = 0
        self.medianQtlLen = 0
        self.longestQtlLen = 0
        self.shortestQtlLen = 0
        self.totalAnnoGenes = 0
        # self.totalTargetGenes = 0
        # self.totalEnrichTargetGenes = 0
        # self.aveEnrichGenes = 0

    def getQtlLengthInfo(self, qtlLenSeries):
        self.medianQtlLen = np.median(qtlLenSeries)
        self.longestQtlLen = np.max(qtlLenSeries)
        self.shortestQtlLen = np.min(qtlLenSeries)

    def getTotalAnnoGenes(self, qtlAllGenesSeries):
        self.totalAnnoGenes = len(self.mergeSeries2List(qtlAllGenesSeries))

    # def getTotalTargetGenes(self, qtlAllMetaGenesSeries):
    #     self.totalTargetGenes = len(self.mergeSeries2List(qtlAllMetaGenesSeries))

    # def getEnrichTargetGenes(self, myDataFrame):
    #     filteredSeries = myDataFrame.loc[myDataFrame.qtl_p_value <= 0.05, "qtl_meta_gene"]
    #     self.totalEnrichTargetGenes = len(self.mergeSeries2List(filteredSeries))

    def getInit(self, tableData):
        self.qtlDetected = len(tableData)
        self.getQtlLengthInfo(tableData.qtl_length)
        self.getTotalAnnoGenes(tableData.qtl_all_gene)
        # self.getTotalTargetGenes(tableData.qtl_meta_gene)
        # self.getEnrichTargetGenes(tableData)

    def mergeSeries2List(self, mySeries, sep=";"):
        return set().union(*mySeries.apply(lambda x: str(x).split(sep)).to_list())


class AllTraitStatistics():
    def __init__(self):
        self.totalTraits = 0
        self.filteredTraits = 0
        self.clusteredTraits = 0
        self.unclusteredTraits = 0

    def getFilteredTraits(self, traitSeries):
        self.filteredTraits = len(self.mergeSeries2List(traitSeries, sep=","))
        self.totalTraits = len(self.mergeSeries2List(traitSeries, sep=","))

    def getClusteredTraits(self, labelSeries):
        counter = Counter(labelSeries.to_list())
        clustered = [i for i in counter if counter[i] != 1]
        unclustered = [i for i in counter if counter[i] == 1]
        self.clusteredTraits = len(clustered)
        self.unclusteredTraits = len(unclustered)

    # def getUnclusterTraits(self, labelSeries):
    #     self.clusteredTraits = len(set(labelSeries.to_list()))

    def getInit(self, tableData):
        self.getFilteredTraits(tableData.phe_name)
        self.getClusteredTraits(tableData.phe_name)
        # self.getUnclusterTraits(tableData.phe_name)

    def mergeSeries2List(self, mySeries, sep=";"):
        return set().union(*mySeries.apply(lambda x: str(x).split(sep)).to_list())


class SingleQtlStatistics():
    def __init__(self, myRow, rowIndex):
        self.qtlName = "QTL{}_chr{}-{}-{}".format(rowIndex, myRow.CHR, myRow.qtl_start, myRow.qtl_end)
        self.qtlPosition = "{}:{}-{}".format(myRow.CHR, myRow.qtl_start, myRow.qtl_end)
        self.peakSNP = "{}, {}".format(myRow.SNP, myRow.SNP)
        self.traitNames = myRow.phe_name.split(",")
        self.totalGenesInQtl = str(myRow.qtl_all_gene).split(";")
        # self.targetGenesInQtl = str(myRow.qtl_meta_gene).split(";")
        # self.enrichTargetGenes = 0
        self.pvalue = myRow.P


################# Functions ##################
def resolveDir(dirName):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    os.chdir(dirName)


def getGeneFunc(myGeneFuncFile, sep=","):
    gene2func = pd.read_csv(myGeneFuncFile, sep=sep)
    tmpDict = gene2func.to_dict("index")
    geneFuncDict = {}
    for i in tmpDict:
        geneFuncDict[tmpDict[i]["geneid"]] = tmpDict[i]
    return geneFuncDict


def getAllQtlSummary(allQtlStatistics, doc=None, tag=None, text=None, line=None):
    with tag("div", id="qtlSummary"):
        line("h1", "Summary of QTLs detected by local GWAS", style="text-align: center;")
        multi_trait_File = "./Manhattan.multi_trait.jpg"
        with tag("div", style="text-align: center;margin-top: 50px;"):
            doc.stag("img", klass="img-fluid", src=multi_trait_File, alt="multi_trait")
        with tag("div", style="margin-top: 100px;"):
            with tag("div"):
                line("h2", "QTL summarization criteria")
                line("p", "The cutoffs in generating and filtering QTLs")
            with tag("div"):
                line("h2", "QTL statistics table")
                with tag("div", klass="table-responsive"):
                    with tag("table", klass="table"):
                        with tag("thead"):
                            with tag("tr", klass="table-success"):
                                line("th", "Categories", style="width: 50%")
                                line("th", "Statistics value", style="width: 50%")
                        with tag("tbody"):
                            with tag("tr"):
                                line("td", "Genotype file")
                                line("td", "XXXXXX")
                            with tag("tr"):
                                line("td", "Number of SNPs for local GWAS")
                                line("td", str(allQtlStatistics.totalSNPs))
                            with tag("tr"):
                                line("td", "Number of detected QTLs")
                                line("td", str(allQtlStatistics.qtlDetected))
                            with tag("tr"):
                                line("td", "Median QTL length")
                                line("td", str(allQtlStatistics.medianQtlLen))
                            with tag("tr"):
                                line("td", "Longest QTL length")
                                line("td", str(allQtlStatistics.longestQtlLen))
                            with tag("tr"):
                                line("td", "Shortest QTL length")
                                line("td", str(allQtlStatistics.shortestQtlLen))
                            with tag("tr"):
                                line("td", "Total genes in the QTLs")
                                line("td", str(allQtlStatistics.totalAnnoGenes))
                            # with tag("tr"):
                            #     line("td", "Total target genes in the QTLs")
                            #     line("td", str(allQtlStatistics.totalTargetGenes))
                            # with tag("tr"):
                            #     line("td", "Total enrichment of target genes")
                            #     line("td", str(allQtlStatistics.totalEnrichTargetGenes))
                            # with tag("tr"):
                            #     line("td", "Average enrichment of target genes")
                            #     line("td", str(allQtlStatistics.aveEnrichGenes))


def getAllTraitSummary(allTraitStatistics, doc=None, tag=None, text=None, line=None):
    # doc, tag, text, line = Doc().ttl()
    with tag("div", id="traitSummary"):
        line("h1", "Summary of omics traits detected by local GWAS", style="text-align: center;")
        # heatmapFile = "/home/xufeng/xufeng/Projects/MODAS/xufeng1/assets/img/heatmap.png"
        # with tag("div", style="text-align: center;margin-top: 50px;"):
        #     doc.stag("img", klass="img-fluid", src=heatmapFile, alt="heatmap")
        with tag("div", style="margin-top: 100px;"):
            with tag("div"):
                line("h2", "Omics trait filtration criteria")
                line("p", "Parameters and pipelines in filtering traits")
            with tag("div"):
                line("h2", "Trait statistics table")
                with tag("div", klass="table-responsive"):
                    with tag("table", klass="table"):
                        with tag("thead"):
                            with tag("tr", klass="table-success"):
                                line("th", "Categories", style="width: 50%")
                                line("th", "Statistics value", style="width: 50%")
                        with tag("tbody"):
                            with tag("tr"):
                                line("td", "Omics trait type")
                                line("td", "Metabolome")
                            with tag("tr"):
                                line("td", "Total number of traits")
                                line("td", str(allTraitStatistics.totalTraits))
                            with tag("tr"):
                                line("td", "Filtered number of traits")
                                line("td", str(allTraitStatistics.filteredTraits))
                            with tag("tr"):
                                line("td", "Number of clustered traits")
                                line("td", str(allTraitStatistics.clusteredTraits))
                            # with tag("tr"):
                            #     line("td", "Number of modules of clustered traits")
                            #     line("td", "Cell 2")
                            with tag("tr"):
                                line("td", "Number of unclustered traits")
                                line("td", str(allTraitStatistics.unclusteredTraits))
    # return doc.getvalue()


def getSingleQtlInfo(singleQtl, index, geneFuncDict, doc=None, tag=None, text=None, line=None):
    with tag("div", klass="container"):
        with tag("div", klass="row"):
            with tag("div", klass="col"):
                line("h1", "Summary in " + singleQtl.qtlName, style="text-align: center;")
                with tag("div", style="margin-top: 100px;"):
                    with tag("div"):
                        line("h2", "QTL summarization criteria for whole -genome GWAS")
                        line("p", "Parameters and pipeline for filtering traits")
                    with tag("div"):
                        line("h2", "QTL statistics table")
                        with tag("div", klass="table-responsive"):
                            with tag("table", klass="table"):
                                with tag("thead"):
                                    with tag("tr", klass="table-success"):
                                        line("th", "QTL", style="width: 50%")
                                        line("th", "Statistics", style="width: 50%")
                                with tag("tbody"):
                                    with tag("tr"):
                                        line("td", "QTL position")
                                        line("td", str(singleQtl.qtlPosition))
                                    with tag("tr"):
                                        line("td", "Peak SNP ID and position")
                                        line("td", str(singleQtl.peakSNP))
                                    with tag("tr"):
                                        line("td", "Number of total genes in the QTL")
                                        line("td", str(len(singleQtl.totalGenesInQtl)))
                                    # with tag("tr"):
                                    #     line("td", "Number of target genes in the QTL")
                                    #     line("td", str(len(singleQtl.targetGenesInQtl)))
                                    # with tag("tr"):
                                    #     line("td", "Enrichment of target genes")
                                    #     line("td", str(singleQtl.enrichTargetGenes))
                                    with tag("tr"):
                                        line("td", "Enrichment significance vs background")
                                        line("td", str(singleQtl.pvalue))
                    # with tag("div"):
                    #     line("h2", "List of target genes in the QTL")
                    #     with tag("div", klass="table-responsive"):
                    #         with tag("table", klass="table"):
                    #             with tag("thead"):
                    #                 with tag("tr", klass="table-success"):
                    #                     line("th", "Gene ID", style="width: 25%")
                    #                     line("th", "Alias ID", style="width: 25%")
                    #                     line("th", "Position", style="width: 25%")
                    #                     line("th", "Function", style="width: 25%")
                    #             with tag("tbody"):
                    #                 for qtl in singleQtl.targetGenesInQtl:
                    #                     if qtl == "nan":
                    #                         continue
                    #                     with tag("tr"):
                    #                         line("td", str(geneFuncDict[qtl]["geneId"].strip()))
                    #                         line("td", str(geneFuncDict[qtl]["aliasId"].strip()))
                    #                         line("td", str(geneFuncDict[qtl]["position"].strip()))
                    #                         line("td", str(geneFuncDict[qtl]["function"].strip()))
                    with tag("div"):
                        line("h2", "List of total genes in the QTL")
                        with tag("div", klass="table-responsive"):
                            with tag("table", klass="table"):
                                with tag("thead"):
                                    with tag("tr", klass="table-success"):
                                        line("th", "Gene ID", style="width: 25%")
                                        line("th", "Alias ID", style="width: 25%")
                                        line("th", "Position", style="width: 25%")
                                        line("th", "Function", style="width: 25%")
                                with tag("tbody"):
                                    for qtl in singleQtl.totalGenesInQtl:
                                        if qtl == "nan":
                                            continue
                                        with tag("tr"):
                                            line("td", str(geneFuncDict[qtl]["geneid"].strip()))
                                            line("td", str(geneFuncDict[qtl]["aliasid"].strip()))
                                            line("td", str(geneFuncDict[qtl]["position"].strip()))
                                            line("td", str(geneFuncDict[qtl]["function"].strip()))


# def getSingleTrait(traitName, doc=None, tag=None, text=None, line=None):
#     with tag("div", klass="container"):
#         with tag("div", klass="row"):
#             with tag("div", klass="col"):
#                 line("h1", "Details in " + traitName, style="text-align: center;")
#                 doc.stag("img", klass="img-fluid", src="../../assets/img/manhattan.jpg")
#                 with tag("div", klass="row"):
#                     with tag("div", klass="col-6"):
#                         doc.stag("img", klass="img-fluid", src="../../assets/img/qqplot.png")
#                     with tag("div", klass="col-6"):
#                         doc.stag("img", klass="img-fluid", src="../../assets/img/boxplot.png")

def getListItem(data, qtlName=None, traitName=None, doc=None, tag=None, text=None, line=None, mainPage=False):
    for index, row in data.iterrows():
        qtlItem = SingleQtlStatistics(row, index)
        if qtlName and qtlName == qtlItem.qtlName:
            expand = "true"
            faPlusOrMinus = "fa-minus"
            myClass = "list-unstyled collapse nav nav-pills show"
            active = " active"
        else:
            expand = "false"
            faPlusOrMinus = "fa-plus"
            myClass = "list-unstyled collapse nav nav-pills"
            active = ""
        if mainPage:
            relativeDir = os.path.join("", qtlItem.qtlName)
        else:
            if qtlName == qtlItem.qtlName:
                relativeDir = ""
            else:
                relativeDir = os.path.join("../", qtlItem.qtlName)
        with tag("li"):
            with tag("div", klass="qtlItem" + active):
                with tag("a", ("href", os.path.join(relativeDir, qtlItem.qtlName + ".html")), klass="qtlLink"):
                    text(qtlItem.qtlName)
                with tag("a", ("href", "#" + qtlItem.qtlName), ("data-toggle", "collapse"), ("aria-expanded", expand)):
                    line("i", "", klass="fa " + faPlusOrMinus)
            with tag("ul", ("class", myClass), ("id", qtlItem.qtlName), ("aria-expanded", expand)):
                for i in qtlItem.traitNames:
                    with tag("li"):
                        href = os.path.join(relativeDir, i + ".html")
                        if traitName and traitName == i:
                            with tag("a", ("href", href), ("class", "active"), ("aria-selected", "true")):
                                line("i", "", klass="fa fa-link")
                                text(" " + i)
                        else:
                            with tag("a", ("href", href), ("aria-selected", "false")):
                                line("i", "", klass="fa fa-link")
                                text(" " + i)


def generateMainPage(data, allQtlStatistics, allTraitStatistics):
    doc, tag, text, line = Doc().ttl()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            doc.stag('meta', charset='utf-8')
            doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
            line('title', 'MODAS main page')
            doc.stag('link', rel='stylesheet', href='assets/bootstrap/css/bootstrap.min.css')
            doc.stag('link', rel='stylesheet', href='assets/fonts/font-awesome.min.css')
            doc.stag('link', rel='stylesheet', href="assets/css/modas.css")
            # doc.stag('link', rel='stylesheet', href="assets/css/styles.css")

        with tag("body"):
            with tag("div", id="sidebar-test"):
                with tag("div", klass="sidebar-header"):
                    with tag("h2"):
                        line("a", "MODAS", href="mainPage.html", klass="modas")
                with tag("ul"):
                    getListItem(data, doc=doc, tag=tag, text=text, line=line, mainPage=True)

            with tag("div", klass="content"):
                with tag("div", klass="container"):
                    with tag("div", klass="row"):
                        with tag("div", klass="col"):
                            getAllQtlSummary(allQtlStatistics, doc=doc, tag=tag, text=text, line=line)
                            doc.stag("hr", style="margin-bottom: 50px;margin-top: 50px;")
                            getAllTraitSummary(allTraitStatistics, doc=doc, tag=tag, text=text, line=line)
            line("script", "", src="assets/js/jquery.min.js")
            line("script", "", src="assets/bootstrap/js/bootstrap.min.js")
            line("script", "", src="assets/js/modas.js")

    mainPageOut = open("mainPage.html", "w")
    res = indent(doc.getvalue(), indentation="    ")
    print(res, file=mainPageOut)
    mainPageOut.close()


def generateSingleQtlPage(data, geneFuncDict):
    for index, row in data.iterrows():
        qtlItem = SingleQtlStatistics(row, index)
        if not os.path.exists(qtlItem.qtlName):
            os.makedirs(qtlItem.qtlName)
        out = open(os.path.join(qtlItem.qtlName, qtlItem.qtlName + ".html"), "w")

        doc, tag, text, line = Doc().ttl()
        doc.asis('<!DOCTYPE html>')
        with tag('html'):
            with tag('head'):
                doc.stag('meta', charset='utf-8')
                doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
                line('title', 'Summary information in QTL ' + qtlItem.qtlName)
                doc.stag('link', rel='stylesheet', href='../assets/bootstrap/css/bootstrap.min.css')
                doc.stag('link', rel='stylesheet', href='../assets/fonts/font-awesome.min.css')
                doc.stag('link', rel='stylesheet', href="../assets/css/modas.css")
                # doc.stag('link', rel='stylesheet', href="assets/css/styles.css")
            with tag("body"):
                with tag("div", id="sidebar-test"):
                    with tag("div", klass="sidebar-header"):
                        with tag("h2"):
                            line("a", "MODAS", href="../mainPage.html", klass="modas")
                    with tag("ul"):
                        getListItem(data, qtlName=qtlItem.qtlName, doc=doc, tag=tag, text=text, line=line)

                with tag("div", klass="content"):
                    getSingleQtlInfo(qtlItem, index, geneFuncDict, doc=doc, tag=tag, text=text, line=line)

                line("script", "", src="../assets/js/jquery.min.js")
                line("script", "", src="../assets/bootstrap/js/bootstrap.min.js")
                line("script", "", src="../assets/js/modas.js")

                customJs = '''
                    <script>
                        var offestFromTop = %d * 45 + 68;
                        $('#sidebar-test').scrollTop(offestFromTop);

                        function clickItem(event) {
                            var target = event.currentTarget;
                            $(target).parent().removeClass(".active").addClass(".active");

                            var index = $("div.qtlItem").index($(this).parent());
                            var offestFromTop = index * 45 + 68;
                            $('#sidebar-test').scrollTop(offestFromTop);
                        }
                        if ($("div.qtlItem .qtlLink")) {
                            var qtlLink = $("div.qtlItem .qtlLink");
                            for (var i = 0; i < qtlLink.length; i++) {
                                var item = qtlLink[i];
                                item.onclick = clickItem;
                            }
                        }
                    </script>
                ''' % (index)
                doc.asis(customJs)

        res = indent(doc.getvalue(), indentation="    ")
        print(res, file=out)
        out.close()


def generateSingleTraitPage(data, manhattanDir, qqDir, boxplotDir):
    for index, row in data.iterrows():
        qtlItem = SingleQtlStatistics(row, index)
        for traitName in qtlItem.traitNames:
            # getListItem(qtlItem)
            out = open(os.path.join(qtlItem.qtlName, traitName + ".html"), "w")
            doc, tag, text, line = Doc().ttl()
            doc.asis('<!DOCTYPE html>')
            with tag('html'):
                with tag('head'):
                    doc.stag('meta', charset='utf-8')
                    doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0, shrink-to-fit=no')
                    line('title', 'Detailed information in trait ' + traitName)
                    doc.stag('link', rel='stylesheet', href='../assets/bootstrap/css/bootstrap.min.css')
                    doc.stag('link', rel='stylesheet', href='../assets/fonts/font-awesome.min.css')
                    doc.stag('link', rel='stylesheet', href="../assets/css/modas.css")
                    # doc.stag('link', rel='stylesheet', href="assets/css/styles.css")
                with tag("body"):
                    with tag("div", id="sidebar-test"):
                        with tag("div", klass="sidebar-header"):
                            with tag("h2"):
                                line("a", "MODAS", href="../mainPage.html", klass="modas")
                        with tag("ul"):
                            getListItem(data, qtlItem.qtlName, traitName, doc=doc, tag=tag, text=text, line=line)

                    manhattanFile = glob.glob(os.path.join(manhattanDir, "Manhattan.*_{}.jpg".format(traitName)))[0]
                    boxplotFile = glob.glob(os.path.join(boxplotDir, "{}_*.jpg".format(traitName)))[0]
                    qqplotFile = glob.glob(os.path.join(qqDir, "QQplot.*_{}.jpg".format(traitName)))[0]
                    with tag("div", klass="content"):
                        with tag("div", klass="container"):
                            with tag("div", klass="row"):
                                with tag("div", klass="col"):
                                    line("h1", "Details in " + traitName, style="text-align: center;")
                                    doc.stag("img", klass="img-fluid", src='../'+'/'.join(manhattanFile.split('/')[-2:]))
                                    with tag("div", klass="row"):
                                        with tag("div", klass="col-6"):
                                            doc.stag("img", klass="img-fluid", src='../'+'/'.join(qqplotFile.split('/')[-2:]))
                                        with tag("div", klass="col-6"):
                                            doc.stag("img", klass="img-fluid", src='../'+'/'.join(boxplotFile.split('/')[-2:]))

                    line("script", "", src="../assets/js/jquery.min.js")
                    line("script", "", src="../assets/bootstrap/js/bootstrap.min.js")
                    line("script", "", src="../assets/js/modas.js")

                    customJs = '''
                        <script>
                            var offestFromTop = %d * 45 + 68;
                            $('#sidebar-test').scrollTop(offestFromTop);

                            function clickItem(event) {
                                var target = event.currentTarget;
                                $(target).parent().removeClass(".active").addClass(".active");

                                var index = $("div.qtlItem").index($(this).parent());
                                var offestFromTop = index * 45 + 68;
                                $('#sidebar-test').scrollTop(offestFromTop);
                            }
                            if ($("div.qtlItem .qtlLink")) {
                                var qtlLink = $("div.qtlItem .qtlLink");
                                for (var i = 0; i < qtlLink.length; i++) {
                                    var item = qtlLink[i];
                                    item.onclick = clickItem;
                                }
                            }
                        </script>
                    ''' % (index)
                    doc.asis(customJs)

            res = indent(doc.getvalue(), indentation="    ")
            print(res, file=out)
            out.close()


def generateHtml(qtl_anno, myFuncFile, out_dir, totalSNPs):
    manhattanDir = os.path.abspath(out_dir+'/manhattan_plot')
    qqplotDir = os.path.abspath(out_dir+'/qqplot')
    boxplotDir = os.path.abspath(out_dir+'/boxplot')

    allQtlStatistics = AllQtlStatistics()
    allQtlStatistics.getInit(qtl_anno)
    allQtlStatistics.totalSNPs = totalSNPs
    allTraitStatistics = AllTraitStatistics()
    allTraitStatistics.getInit(qtl_anno)

    geneFuncDict = getGeneFunc(myFuncFile, "\t")

    resolveDir(out_dir)
    generateMainPage(qtl_anno, allQtlStatistics, allTraitStatistics)
    generateSingleQtlPage(qtl_anno, geneFuncDict)

    generateSingleTraitPage(qtl_anno, manhattanDir, qqplotDir, boxplotDir)

