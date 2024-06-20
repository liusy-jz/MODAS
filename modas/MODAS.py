#!/usr/bin/env python
import sys, time
import pandas as pd
import numpy as np
import datetime
import argparse
import traceback
import subprocess
import warnings
import shutil
import glob
import os
import re


warnings.filterwarnings("ignore")


class Logger(object):

    def __init__(self, fn):
        self.log_file = open(fn, 'w')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.

        '''
        if isinstance(msg,list):
            msg = ' '.join(msg)
        else:
            msg = ('{T}: '+msg).format(T=time.ctime())
        print(msg,file=self.log_file)
        print(msg)


class ArgsError(Exception):
    def myexcepthook(type, value, tb):
        msg = ''.join(traceback.format_exception(type, value, tb))
        print(msg, end = '')

    sys.excepthook = myexcepthook


class FileTypeError(Exception):
    def myexcepthook(type, value, tb):
        msg = ''.join(traceback.format_exception(type, value, tb))
        print(msg, end = '')

    sys.excepthook = myexcepthook


def genoidx(args, log):
    import modas.genoidx as gi

    if not os.path.exists(args.g) and not os.path.exists(args.g+'.bed'):
        raise ArgsError('the genotype file {} is not exists.'.format(args.g))
    if args.convert:
        log.log('begin convert hapmap genotype file to plink bed format genotype file.')
        b = gi.hap2plink_bed(args.g, args.o)
        if b:
            args.g = args.o
            log.log('Successfully convert hapmap genotype file to plink bed format genotype file.')
        else:
            raise FileTypeError('{} is not a hapmap file or {} is empty.'.format(args.g, args.g))
    g = gi.read_genotype(args.g)
    if g is None:
            raise FileTypeError('{0}.bed is not a plink bed format file or {0}.bim file and {0}.fam file is not in the same folder with {0}.bed, or {0}.bed is empty.'.format(args.g))
    if args.clump:
        log.log('begin snp clumping...')
        if os.path.exists(args.o + '_clump.bed'):
            os.remove(args.o + '_clump.bed')
            os.remove(args.o + '_clump.bim')
            os.remove(args.o + '_clump.fam')
        g_clump_list = gi.snp_clumping(args.g+'.bed', 0.8, args.o + '_clump.bed')
        log.log('snp clumping done.')
        log.log('There are {0} snps in input genotype file, after clumping, {1} snp left.'.format(g_clump_list[0],g_clump_list[1]))
        log.log('clumped genotype file is saved as {0}_clump.bed, {0}_clump.bim and {0}_clump.fam.'.format(args.o))
    if args.genome_cluster:
        if args.w <= 0:
            raise ArgsError('the window size of chromosome should larger than zero')
        if args.s < 0:
            raise ArgsError('the step size should larger than zero')
        log.log('begin genome cluster analysis...')
        if args.clump:
            g_clump = gi.read_genotype(args.o + '_clump')
        else:
            g_clump = g
        #genome_cluster,variant_cluster = gi.genome_cluster(g_clump,args.w,args.s,args.p)
        genome_cluster = gi.genome_cluster(g_clump, args.w, args.s, args.p)
        if genome_cluster.empty:
            log.log('genome cluster analysis failed,Please check genotype file.')
            sys.exit()
        genome_cluster.to_csv(args.o+'.genome_cluster.csv')
        # variant_cluster.to_csv(args.o+'.variant_cluster.csv',index=False)
    log.log('genoidx is finishd!')


def phenorm(args, log):
    import modas.genoidx as gi
    import modas.phenorm as pn

    log.log('begin phenotype filter/normalize analysis...')
    if not os.path.exists(args.phe):
        raise ArgsError('the phenotype file {} is not exists.'.format(args.phe))
    phe = pd.read_csv(args.phe, index_col=0)
    if phe.empty:
        log.log('the phenotype file {} is empty and software is terminated.'.format(args.phe))
        sys.exit()
    if args.v is None and args.r is None and not args.log2 and not args.log10 and not args.ln and not args.norm and not args.qqnorm:
        log.log('No phenotype filter process or normalize process choosed, please select one of the argument in -a, -r, -log2, -log10, -ln, -norm, -qqnorm.')
        sys.exit()
    log.log('Read phenotype file finished,There are {} phenotypes in phenotype file.'.format(phe.shape[1]))
    if args.r is not None:
        if args.r < 0 or args.r > 1:
            raise ArgsError('the interval of missing ration is between 0 and 1.')
        log.log('filter phenotype missing ratio larger than {}'.format(args.r))
        phe = pn.missing_filter(phe, args.r)
        log.log('after filtered, There are {} phenotypes left.'.format(phe.shape[1]))
    if args.v is not None:
        if args.v <= 0:
            raise ArgsError('metabolite abundance cutoff should larger than zero.')
        log.log('filter phenotype values less than {}'.format(args.v))
        phe = pn.abundance_filter(phe, args.v)
        log.log('after filtered, There are {} phenotypes left.'.format(phe.shape[1]))
    phe = phe.loc[:, (phe.var()-0) >= 1e-6]
    if args.pca:
        log.log('correct trait by PCA')
        g = gi.read_genotype(args.g)
        if g is None:
            raise FileTypeError('{0}.bed is not a plink bed format file or {0}.bim file and {0}.fam file is not in the same folder with {0}.bed, or {0}.bed is empty.'.format(args.g))
        pca = pn.pc_calc(args.g + '.bed', 3)
        if pca is None:
            log.log('calculate PCA failed, phenotype correction analysis is not performed')
        else:
            pca_phe_idx = pca.index.intersection(phe.index)
            phe = pn.trait_correct(pca.reindex(pca_phe_idx), phe.reindex(pca_phe_idx))
            log.log('correcting phenotype is successed')
    if args.log2:
        log.log('log2 scale phenotype values.')
        phe = pn.log2_scale(phe)
    if args.log10:
        log.log('log10 scale phenotype values.')
        phe = pn.log10_scale(phe)
    if args.ln:
        log.log('ln scale phenotype values.')
        phe = pn.ln_scale(phe)
    if args.norm:
        log.log('normalize phenotype values')
        phe = pn.normalize_scale(phe)
    if args.qqnorm:
        log.log('normal quantile transform phenotype values')
        phe = phe.apply(pn.qqnorm)
    phe.to_csv(args.o+'.normalized_phe.csv', na_rep='NA')
    log.log('phenotype filter/normalize analysis is finished.')


def prescreen(args, log):
    import modas.prescreen as ps

    log.log('begin phenotype prescreen analysis...')
    if not os.path.exists(args.phe):
        raise ArgsError('the phenotype file {} is not exists.'.format(args.phe))
    phe = pd.read_csv(args.phe, index_col=0)
    if phe.empty:
        log.log('the phenotype file {} is empty and software is terminated.'.format(args.phe))
        sys.exit()
    log.log('Read phenotype file finished,There are {} phenotypes in phenotype file.'.format(phe.shape[1]))
    if not os.path.exists(args.genome_cluster):
        raise ArgsError('genome cluster file is not exists, can not run phenotype prescreen analysis,please assign a exist file for -genome_cluster.')
    if args.gwas_model == 'MLM':
        if args.g is None:
            log.log('the plink bed format genotype  not provided, Using LM(linear model) for pseudo-genotype file GWAS analysis.')
            args.gwas_model = 'LM'
        elif not os.path.exists(args.g+'.bed'):
            log.log('the plink bed format genotype file {}.bed is not exists. Using LM(linear model) for pseudo-genotype file GWAS analysis'.format(args.g))
            args.gwas_model = 'LM'
    genome_cluster = pd.read_csv(args.genome_cluster, index_col=0)
    if genome_cluster.empty:
        log.log('genome cluster file {} is empty and software is terminated.'.format(args.genome_cluster))
        sys.exit()
    if os.path.exists('tmp_omics_phe_bimbam'):
        shutil.rmtree('tmp_omics_phe_bimbam')
    os.mkdir('tmp_omics_phe_bimbam')
    if os.path.exists('output'):
        os.rename('output','output'+datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))
    #if args.gwas_model =='MLM':
    #    fam = pd.read_csv(args.g+'.fam', sep=r'\s+', header=None,index_col=0)
    #    sample_id = fam.index.intersection(genome_cluster.index)
    #    genome_cluster = genome_cluster.reindex(sample_id)
    #    geno_prefix = args.g.split('/')[-1]
    #else:
    #    sample_id = genome_cluster.index
    #    geno_prefix = '.'.join(args.genome_cluster.split('/')[-1].split('.')[:-2])
    #phe = phe.reindex(sample_id)
    phe = phe.reindex(genome_cluster.index)
    geno_prefix = '.'.join(args.genome_cluster.split('/')[-1].split('.')[:-2])
    log.log('convert genome cluster file to bimbam genotype file.')
    a, g = ps.qtl_pc2bimbam(genome_cluster)
    a.to_csv('tmp_omics_phe_bimbam/'+geno_prefix+'_qtl_pc.anno.txt',index=False,header=None)
    g.to_csv('tmp_omics_phe_bimbam/'+geno_prefix+'_qtl_pc.geno.txt',index=False,header=None)
    if args.gwas_suggest_pvalue is None:
        args.gwas_suggest_pvalue = 1.0 / genome_cluster.shape[1]
    log.log('running pseudo-genotype GWAS for significant phenotype filter.')
    s = ps.qtl_pc_gwas_parallel(phe, 'tmp_omics_phe_bimbam', args.p, args.g, geno_prefix, args.gwas_model)
    # if args.lm_suggest_pvalue<1:
    #     log.log('run lm model for genome cluster pseudo-SNP filter.')
    #     s = ps.qtl_pc_lm_gwas_parallel(phe,'tmp_omics_phe_bimbam',args.p,args.g)
    #     ps.generate_omics_qtl_pc_bimbam(phe,a,g,args.lm_suggest_pvalue,args.p)
    #     log.log('run lmm model for Significant phenotype filter.')
    #     s = ps.qtl_pc_lmm_gwas_parallel(phe,'tmp_omics_phe_bimbam',args.p,args.g,sample_id)
    # else:
    #     log.log('run lmm model for Significant phenotype filter.')
    #     s = ps.qtl_pc_lmm_gwas_parallel(phe, 'tmp_omics_phe_bimbam', args.p, args.g, sample_id, args.lm_suggest_pvalue)
    sig_omics_phe, phe_sig_qtl = ps.prescreen(phe, args.gwas_suggest_pvalue)
    sig_omics_phe.to_csv(args.o+'.sig_omics_phe.csv', na_rep='NA')
    phe_sig_qtl.to_csv(args.o+'.phe_sig_qtl.csv', index=False)
    log.log('prescreened phenotype saved successfully')


def regiongwas(args, log):
    import modas.regiongwas as rg

    log.log('begin regiongwas analysis...')
    if not os.path.exists(args.phe):
        raise ArgsError('the phenotype file {} is not exists.'.format(args.phe))
    phe = pd.read_csv(args.phe, index_col=0)
    if phe.empty:
        log.log('the phenotype file {} is empty and software is terminated.'.format(args.phe))
        sys.exit()
    log.log('Read phenotype file finished,There are {} phenotypes in phenotype file.'.format(phe.shape[1]))
    if not os.path.exists(args.phe_sig_qtl):
        raise ArgsError('the phenotype significant qtl file {} is not exists.'.format(args.phe_sig_qtl))
    phe_sig_qtl = pd.read_csv(args.phe_sig_qtl)
    if phe_sig_qtl.empty:
        log.log('the phenotype significant qtl file {} is empty and software is terminated.'.format(args.phe))
        sys.exit()
    if not os.path.exists(args.g+'.bed'):
        raise ArgsError('the plink bed format genotype file {} is not exists.'.format(args.g))
    if os.path.exists('tmp_qtl_bed'):
        shutil.rmtree('tmp_qtl_bed')
    os.mkdir('tmp_qtl_bed')
    if os.path.exists('tmp_rs_dir'):
        shutil.rmtree('tmp_rs_dir')
    os.mkdir('tmp_rs_dir')
    log.log('begin generate phenotype significant gwas region SNP genotype file...')
    rg.generate_qtl_batch(phe,phe_sig_qtl,args.g,args.p,'tmp_qtl_bed','tmp_rs_dir')
    log.log('begin region gwas analysis of phenotype significant gwas region...')
    s = rg.region_gwas_parallel('tmp_qtl_bed', args.p,args.g,args.gwas_model)
    log.log('begin generate QTL of phenotypes...')
    s = rg.generate_clump_input('output', args.gwas_model)
    if args.p1 is None:
        bim = pd.read_csv(args.g+'.bim', header=None, sep='\s+')
        args.p1 = 1.0 / bim.shape[0]
    if args.p2 is None:
        args.p2 = args.p1 * 10
    s = rg.plink_clump('tmp_qtl_bed', args.p1, args.p2, args.p)
    qtl_res, bad_qtl = rg.generate_qtl('clump_result', args.p2)
    qtl_res.index = np.arange(qtl_res.shape[0])
    qtl_res.to_csv(args.o+'.region_gwas_qtl_res.csv', index=False)
    bad_qtl.to_csv(args.o+'.region_gwas_bad_qtl_res.csv', index=False)
    log.log('generate QTL done.')
    if args.cluster:
        log.log('begin cluster phenotype analysis...')
        clustered_phe, phe_labeled = rg.phe_PCA(phe, qtl_res)
        clustered_phe.to_csv(args.o+'.clustered_phe.csv')
        phe_labeled.to_csv(args.o+'.phe_labeled.csv')
        log.log('cluster phenotype analysis done')
    log.log('region gwas analysis is finished.')


def mr(args, log):
    import modas.genoidx as gi
    import modas.mr as MR

    utils_path = subprocess.check_output('locate modas/utils', shell=True, text=True, encoding='utf-8')
    utils_path = '/'.join(re.search('\n(.*site-packages.*)\n', utils_path).group(1).split('/')[:-1])

    log.log('begin Mendelian Randomization analysis...')
    if not os.path.exists(args.exposure):
        raise ArgsError('the exposure file {} is not exists.'.format(args.exposure))
    mTrait = pd.read_csv(args.exposure, index_col=0)
    mTrait.columns = [col.replace('m/z', 'm.z') for col in mTrait.columns]
    if mTrait.empty:
        log.log('the exposure file {} is empty and software is terminated.'.format(args.exposure))
        sys.exit()
    if not os.path.exists(args.outcome):
        raise ArgsError('the outcome file {} is not exists.'.format(args.outcome))
    pTrait = pd.read_csv(args.outcome, index_col=0)
    if pTrait.empty:
        log.log('the outcome file {} is empty and software is terminated.'.format(args.outcome))
        sys.exit()
    if not os.path.exists(args.qtl):
        raise ArgsError('the exposure QTL file {} is not exists.'.format(args.qtl))
    qtl = pd.read_csv(args.qtl)
    if qtl.empty:
        log.log('the exposure QTL file {} is empty and software is terminated.'.format(args.qtl))
        sys.exit()
    if not os.path.exists(args.g+'.bed'):
        raise ArgsError('the plink bed format genotype file {} is not exists.'.format(args.g))
    if args.lm and args.mlm:
        log.log('please select linear model or mixed linear model for Mendelian Randomization analysis.')
        log.log('software is terminated.')
        sys.exit()
    if not args.lm and not args.mlm:
        log.log('please select a model from linear model and mixed linear model for Mendelian Randomization analysis.')
        log.log('software is terminated.')
        sys.exit()
    if mTrait.columns.isin(qtl.phe_name).sum()>0:
        mTrait = mTrait.loc[:, mTrait.columns.isin(qtl.phe_name)]
    else:
        log.log('there is no mTrait in exposure QTL result, please check your QTL file or exposure file.')
        log.log('software is terminated.')
        sys.exit()
    if args.lm:
        log.log('perform Mendelian Randomization through linear model')
        g = gi.read_genotype(args.g)
        if g is None:
            raise FileTypeError('{0}.bed is not a plink bed format file or {0}.bim file and {0}.fam file is not in the same folder with {0}.bed, or {0}.bed is empty.'.format(args.g))
        g = g.where(g.snp.isin(qtl.SNP), drop=True)
        g = pd.DataFrame(g.values, index=g.sample, columns=g.snp.values)
        ril = g.index.intersection(mTrait.index.intersection(pTrait.index))
        g = g.reindex(ril)
        mTrait = mTrait.reindex(ril)
        pTrait = pTrait.reindex(ril)
        res = MR.MR_parallel(qtl, mTrait, pTrait, g, args.p, args.pvalue)
        res.to_csv(args.o + '.MR.csv', index=False)
        log.log('Successfully perform Mendelian randomization analysis using linear model')
    if args.mlm:
        log.log('perform Mendelian Randomization through mixed linear model')
        MR.generate_geno_batch(qtl, mTrait, pTrait, args.g, args.p, 'tmp_mr_bed', 'tmp_mr_rs')
        MR.calc_MLM_effect('tmp_mr_bed', pTrait, args.p, args.g)
        mTrait_effect, pTrait_effect, pTrait_se = MR.get_MLM_effect_parallell('./output', mTrait, pTrait, args.p)
        res = MR.MR_MLM_parallel(qtl, mTrait_effect, pTrait_effect, pTrait_se, args.p, args.pvalue)
        res.to_csv(args.o+'.MR.csv', index=False)
        log.log('Successfully perform Mendelian randomization analysis using mixed linear model')
    log.log('Mendelian Randomization analysis is successed')
    # res = pd.read_csv('chr_HAMP_kernel_eqtl_new.local_gwas_qtl_res.MR.csv')
    if args.net:
        log.log('Mendelian Randomization network analysis')
        edge_weight = MR.edge_weight(qtl, res)
        edge_weight = edge_weight.loc[edge_weight.weight >= 0.2, :]
        edge_weight = edge_weight.sort_values(by=['mTrait', 'pTrait', 'weight'], ascending=False).drop_duplicates(subset=['mTrait', 'pTrait'])
        edge_weight.to_csv(args.o+'.edgelist', sep='\t', header=None, index=False)
        subprocess.call('java -jar ' + utils_path +'/cluster_one-1.0.jar -f edge_list -F csv '+args.o+'.edgelist'+' >'+args.o+'.cluster_one.result.csv 2>/dev/null',shell=True)
        cluster_one_res = pd.read_csv(args.o+'.cluster_one.result.csv')
        cluster_one_res.loc[(cluster_one_res['P-value'] <= 0.05) & (cluster_one_res['Size'] >= 5), :].to_csv(args.o+'.sig.cluster_one.result.csv', index=False)
        log.log('Mendelian Randomization network analysis is successed')


def visual(args, log):
    import modas.genoidx as gi
    import modas.visual as vis

    utils_path = subprocess.check_output('locate modas/utils', shell=True, text=True, encoding='utf-8')
    utils_path = '/'.join(re.search('\n(.*site-packages.*)\n', utils_path).group(1).split('/')[:-1])

    log.log('begin genome-wide association analysis and visualization')
    if not os.path.exists(args.phe):
        raise ArgsError('the phenotype file {} is not exists.'.format(args.phe))
    phe = pd.read_csv(args.phe, index_col=0)
    if phe.empty:
        log.log('the phenotype file {} is empty and software is terminated.'.format(args.phe))
        sys.exit()
    if not os.path.exists(args.qtl):
        raise ArgsError('the phenotype QTL file {} is not exists.'.format(args.qtl))
    qtl = pd.read_csv(args.qtl)
    if qtl.empty:
        log.log('the phenotype QTL file {} is empty and software is terminated.'.format(args.qtl))
        sys.exit()
    phe = phe.loc[:, phe.columns.isin(qtl.phe_name)]
    if not os.path.exists(args.g + '.bed'):
        raise ArgsError('the plink bed format genotype file {} is not exists.'.format(args.g))
    g = gi.read_genotype(args.g)
    if g is None:
        raise FileTypeError('{0}.bed is not a plink bed format file or {0}.bim file and {0}.fam file is not in the same folder with {0}.bed, or {0}.bed is empty.'.format(args.g))
    log.log('begin GWAS analysis')
    s = vis.gwas(phe, args.g, args.p, args.phe, args.gwas_model)
    if s is None:
        log.log('calculated related_martix faild, please check your genotype file.')
        sys.exit()
    else:
        if not isinstance(s, list) and s == 1:
            log.log('GWAS analysis is failed, please check your genotype file or phenotype file.')
            sys.exit()
        status = np.array(s) == 0
        fail_phe = phe.columns[~status].astype(str)
        if not status.all():
            log.log('There are {0} omics traits GWAS is failed , and omics trait names are {1}'.format(np.sum(~status), ','.join(fail_phe)))
        else:
            log.log('All phenotype GWAS done.')
    if args.visual:
        if not os.path.exists(args.phe):
            raise ArgsError('the gene annotation file {} is not exists.'.format(args.anno))
        anno = pd.read_csv(args.anno, sep='\t')
        if anno.empty:
            log.log('the gene annotation file {} is empty and software is terminated.'.format(args.phe))
            sys.exit()
        qtl_anno = vis.qtl_anno(qtl, anno)
        qtl_anno.to_csv('test_anno.csv', index=False)
        if os.path.exists(args.o):
            shutil.rmtree(args.o)
        os.makedirs(args.o+'/manhattan_plot')
        os.makedirs(args.o+'/qqplot')
        os.makedirs(args.o+'/boxplot')
        shutil.copytree(utils_path + '/assets', args.o + '/assets')
        log.log('begin plot manhattan plot and qqplot')
        s = vis.gwas_plot_parallel(phe, 1e-2, args.p, 'jpg', args.phe, args.gwas_model)
        log.log('begin plot boxplot')
        vis.boxplot(phe, g, qtl)
        log.log('begin plot multi trait manhattan plot')
        vis.multi_trait_plot(phe, './output/', qtl, args.phe, 'multi_trait', 'jpg', args.gwas_model)
        for fn in glob.glob('Manhattan.*jpg'):
            if fn != 'Manhattan.multi_trait.jpg':
                shutil.move(fn, args.o+'/manhattan_plot')
        for fn in glob.glob('QQplot.*jpg'):
            shutil.move(fn, args.o+'/qqplot')
        for fn in glob.glob('*boxplot.jpg'):
            shutil.move(fn, args.o+'/boxplot')
        shutil.move('./Manhattan.multi_trait.jpg', args.o)
        fam = pd.read_csv(args.g+'.fam',sep='\s+')
        vis.generateHtml(qtl_anno, args.anno, args.o, fam.shape[0])
    log.log('Genome-wide association analysis and visualization is down')


def contrast(args, log):
    import modas.genoidx as gi
    import modas.phenorm as pn
    import modas.contrast as ct

    log.log('begin to construct stress response index and identify stress-responsive molecular QTL...')
    if not os.path.exists(args.stress_phe):
        raise ArgsError('the stress omics data file {} is not exists.'.format(args.stress_phe))
    stress_omics = pd.read_csv(args.stress_phe, index_col=0)
    stress_omics.columns = [col.replace('m/z', 'm.z') for col in stress_omics.columns]
    if stress_omics.empty:
        log.log('the stress omics data file {} is empty and software is terminated.'.format(args.stress_phe))
        sys.exit()
    if not os.path.exists(args.control_phe):
        raise ArgsError('the control omics data file {} is not exists.'.format(args.control_phe))
    control_omics = pd.read_csv(args.control_phe, index_col=0)
    control_omics.columns = [col.replace('m/z', 'm.z') for col in control_omics.columns]
    if control_omics.empty:
        log.log('the control omics data file {} is empty and software is terminated.'.format(args.control_phe))
        sys.exit()
    omics_scpca_phe, omics_scpca_test = ct.scpca_omics(stress_omics, control_omics, args.alpha, args.n_comp, args.beta_test_pvalue)
    omics_scpca_phe.to_csv(args.o+'.scpca_pc.phe.csv', na_rep='NA')
    omics_scpca_test.to_csv(args.o+'.scpca_pc.beta_test.csv', index=False)
    log.log('The construction of stress response index is completed.')
    if args.gwas:
        log.log('normal quantile transform phenotype values.')
        omics_scpca_phe = omics_scpca_phe.apply(pn.qqnorm)
        omics_scpca_phe.to_csv(args.o+'.scpca_pc.normalized_phe.csv', na_rep='NA')
        ps_args = argparse.Namespace()
        for arg, value in zip(('phe', 'genome_cluster', 'g', 'gwas_model', 'gwas_suggest_pvalue', 'p','o'), (args.o+'.scpca_pc.normalized_phe.csv', args.genome_cluster, args.g, args.gwas_model, None, args.p, args.o)):
            setattr(ps_args, arg, value)
        rg_args = argparse.Namespace()
        for arg, value in zip(('phe', 'phe_sig_qtl', 'g', 'gwas_model', 'p1', 'p2', 'p', 'cluster', 'o'),(args.o+'.sig_omics_phe.csv', args.o+'.phe_sig_qtl.csv', args.g, args.gwas_model, None, None, args.p, False, args.o)):
            setattr(rg_args, arg, value)
        prescreen(ps_args, log)
        regiongwas(rg_args, log)
        scpca_qtl_res = pd.read_csv(args.o+'.region_gwas_qtl_res.csv')
        geno = gi.read_genotype(args.g)
        var = pd.Series(geno.variant, index=geno.snp)
        geno = geno.sel(variant=var.reindex(scpca_qtl_res.SNP.unique()).values, drop=True)
        geno = pd.DataFrame(geno.values, index=geno.fid, columns=geno.snp)
        scpca_qtl_res = ct.anova(stress_omics, control_omics, scpca_qtl_res, geno)
        scpca_qtl_res.to_csv(args.o+'.region_gwas_qtl_res.anova.csv', index=False)
        os.remove(args.o + '.phe_sig_qtl.csv')
        os.remove(args.o + '.sig_omics_phe.csv')
        os.remove(args.o + '.region_gwas_qtl_res.csv')
        scpca_qtl_res.loc[scpca_qtl_res.pvalue<=0.05, :].to_csv(args.o+'.region_gwas_qtl_res.anova.sig.csv', index=False)
        log.log('two-way ANOVA analysis done')
        log.log('Identification of stress-responsive molecule QTL done.')


def coloc(args, log):
    import modas.genoidx as gi
    import modas.coloc as cl

    log.log('Begin co-localization analysis of multiple traits...')
    if not os.path.exists(args.qtl):
        raise ArgsError('the molecular QTL file {} is not exists.'.format(args.qtl))
    qtl = pd.read_csv(args.qtl)
    if not os.path.exists(args.g+'.bed'):
        raise ArgsError('the genotype file {} is not exists.'.format(args.g))
    g = gi.read_genotype(args.g)
    if not os.path.exists(args.gwas_dir):
        raise ArgsError('the gwas result directory {} is not exists.'.format(args.gwas_dir))
    qtl = cl.qtl_cluster(qtl)
    log.log('Start calculating the kinship matrix for each molecular feature based on the significantly associated SNP')
    kin_info = cl.get_kin_info(qtl, args.gwas_dir, g, args.pvalue)
    if not kin_info:
        log.log('There are no GWAS results for molecular features in the directory or the GWAS results for molecular features do not exceed the threshold')
        log.log('co-localization analysis is terminated.')
        sys.exit()
    log.log('start multi-trait co-localization analysis based on approximate image matching algorithm')
    coloc_res, coloc_pairwise_res, dis = cl.trait_coloc(kin_info, qtl, args.metric, args.eps, args.p)
    coloc_res.to_csv(args.o + '.coloc_res.csv', index=False)
    coloc_pairwise_res.to_csv(args.o + '.coloc_pairwise.csv', index=False)
    dis.to_csv(args.o + '.dis_res.csv', na_rep=1)
    log.log('multi-trait co-localization analysis is down')


def MODAS(args):
    if args.command == 'genoidx':
        log = Logger(args.o+'.genoidx.log.txt')
        log.log(sys.argv)
        genoidx(args, log)
    if args.command == 'phenorm':
        log = Logger(args.o+'.phenorm.log.txt')
        log.log(sys.argv)
        phenorm(args, log)
    if args.command == 'prescreen':
        log = Logger(args.o+'.prescreen.log.txt')
        log.log(sys.argv)
        prescreen(args, log)
    if args.command == 'regiongwas':
        log = Logger(args.o+'.regiongwas.log.txt')
        log.log(sys.argv)
        regiongwas(args, log)
    if args.command == 'mr':
        log = Logger(args.o+'.mr.log.txt')
        log.log(sys.argv)
        mr(args, log)
    if args.command == 'visual':
        log = Logger('.'.join(args.phe.split('/')[-1].split('.')[:-1])+'.visual.log.txt')
        log.log(sys.argv)
        visual(args, log)
    if args.command == 'contrast':
        log = Logger(args.o+'.contrast.log.txt')
        log.log(sys.argv)
        contrast(args, log)
    if args.command == 'coloc':
        log = Logger(args.o+'.coloc.log.txt')
        log.log(sys.argv)
        coloc(args, log)


def main():
    if len(sys.argv) <= 2:
        sys.argv.append('-h')
    args = parser.parse_args()
    MODAS(args)


USAGE = ' MODAS: Multi-omics data association study '

parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)
parser_genoidx = subparsers.add_parser('genoidx', help='generate a pseudo-genotype file', usage='%(prog)s [options]')
parser_genoidx.add_argument('-g', metavar='', help='Genotype file used for genome-wide genotype dimensionality reduction analysis, support plink bed format and hapmap format，the plink bed format file only needs to provide the file prefix, the hapmap format needs to be used together with the parameter -convert')
parser_genoidx.add_argument('-w', default=1000000, type=float, metavar='[default:1000000]', help='Window size for genome-wide genotype dimensionality reduction analysis(bp), default 1000000')
parser_genoidx.add_argument('-s', default=500000, type=float, metavar='[default:500000]', help='Step size for genome-wide genotype dimensionality reduction analysis(bp), default 500000')
parser_genoidx.add_argument('-genome_cluster', action='store_true', help='Perform genome-wide genotype dimensionality reduction analysis on genotype file')
parser_genoidx.add_argument('-convert', action='store_true', help='Convert hapmap genotype format to plink bed format genotype')
parser_genoidx.add_argument('-clump', action='store_true', help='Clumping analysis for genotype file, used to keep a subset of SNPs that are nearly uncorrelated with each other, thereby reducing the number of total SNPs to speed up genome-wide genotype dimensionality reduction analysis')
parser_genoidx.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for genome cluster analysis,max threads is total number of genome chromosomes')
parser_genoidx.add_argument('-o', default='MODAS_genoidx_output', metavar='[default:MODAS_genoidx_output]', help='The prefix of the output file, default MODAS_genoidx_output')
parser_phenorm = subparsers.add_parser('phenorm', help='phenotype normalization and transformation', usage='%(prog)s [options]')
parser_phenorm.add_argument('-phe', metavar='', help='Phenotype file for phenotype preprocessing and transformation, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
parser_phenorm.add_argument('-g', metavar='', help='Genotype file in plink bed format for principal component analysis, used with the parameter -pca to correct differences in phenotype caused by population structure')
parser_phenorm.add_argument('-r', default=None, type=float, metavar='', help='Phenotype missing ratio filter cutoff')
parser_phenorm.add_argument('-v', default=None, type=float, metavar='', help='Phenotype value filter cutoff')
parser_phenorm.add_argument('-log2', action='store_true', help=' log2 transformation of phenotype')
parser_phenorm.add_argument('-log10', action='store_true', help='log10 trasform phenotype')
parser_phenorm.add_argument('-ln', action='store_true', help='ln transform phenotype')
parser_phenorm.add_argument('-norm', action='store_true', help='Boxcox transformation of phenotype')
parser_phenorm.add_argument('-qqnorm', action='store_true', help='Normal quantile transformation of phenotype')
parser_phenorm.add_argument('-pca', action='store_true', help='Correct the differences in phenotype caused by population structure through PCA')
parser_phenorm.add_argument('-o', default='MODAS_phenorm_output', metavar='[default:MODAS_phenorm_output]', help='output file prefix')
parser_prescreen = subparsers.add_parser('prescreen', help='prescreening of associated QTLs with lm/lmm', usage='%(prog)s [options]')
parser_prescreen.add_argument('-g', default=None, metavar='', help='Genotype file in plink bed format, used to calculate the kinship matrix for GWAS analysis')
parser_prescreen.add_argument('-genome_cluster', metavar='', help='Pseudo-genotype file for phenotype pre-screening, generated by subcommand genoidx')
parser_prescreen.add_argument('-phe', metavar='', help='Phenotype file for phenotype pre-screening analysis, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
parser_prescreen.add_argument('-gwas_suggest_pvalue', default=None, type=float, metavar='[default:1/ genome cluster file converted pseudo-genotype number]', help='suggested GWAS P value for pseudo-genotype filtering in regional association analysis, default 1/ pseudo-genotype file pseudo-snp number')
parser_prescreen.add_argument('-gwas_model', default='MLM', type=str, metavar='[default:MLM]', help='model for pseudo-genotype GWAS analysis, supporting LM(linear model), GLM(general linear model) and MLM(mixed linear model) model, default MLM')
parser_prescreen.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for analysis in prescreen sub command')
parser_prescreen.add_argument('-o', default='MODAS_prescreen_output', metavar='[default:MODAS_prescreen_output]', help='The prefix of the output file, default MODAS_prescreen_output')
parser_regiongwas = subparsers.add_parser('regiongwas', help='perform regiongwas to identify QTLs', usage='%(prog)s [options]')
parser_regiongwas.add_argument('-g', metavar='', help='Genotype file in plink bed format, used to calculate the kinship matrix for regional gwas analysis and extract snp in significant association regions of phenotypes')
parser_regiongwas.add_argument('-phe', metavar='', help='Candidate phenotype file generated by subcommand prescreen')
parser_regiongwas.add_argument('-phe_sig_qtl', metavar='', help='Significant association regions of candidate phenotypes file generated by subcommand prescreen')
parser_regiongwas.add_argument('-p1', default=None, type=float, metavar='[default:1/snp number in genotype file]', help='Significance threshold for index SNPs, used to determine the phenotype with QTL and index snps of phenotype')
parser_regiongwas.add_argument('-p2', default=None, type=float, metavar='[default:10/snp number in genotype file]', help='Secondary significance threshold for clumped SNPs, Used to obtain the secondary significant SNP linked to index snp to determine the reliability of the QTL. MODAS outputs the QTL containing 10 secondary significant SNPs as the phenotypic QTL result')
parser_regiongwas.add_argument('-gwas_model', default='MLM', type=str, metavar='[default:MLM]', help='GWAS model for region association analysis, supporting LM, GLM and MLM, default MLM.')
parser_regiongwas.add_argument('-cluster', action='store_true', help='Cluster phenotype by QTL positon')
parser_regiongwas.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for regional gwas analysis, default 1')
parser_regiongwas.add_argument('-o', default='MODAS_regiongwas_out', metavar='[default:MODAS_regiongwas_out]', help='The prefix of the output file, default MODAS_regiongwas_out')
parser_mr = subparsers.add_parser('mr', help='perform Mendelian Randomization and identify function modules', usage='%(prog)s [options]')
parser_mr.add_argument('-g', metavar='', help='Genotype file in plink bed format')
parser_mr.add_argument('-exposure', metavar='', help='Exposure phenotype file, such as mTrait phenotype, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
parser_mr.add_argument('-outcome', metavar='', help='Outcome phenotype file, such as pTrait phenotype, the file format is csv format, the first column and the first row are the names of inbred lines and phenotypes respectively')
parser_mr.add_argument('-qtl', metavar='', help='Exposure phenotype QTL file, generated by subcommand regiongwas')
parser_mr.add_argument('-lm', action='store_true', help='Perform Mendelian Randomization through linear model')
parser_mr.add_argument('-mlm', action='store_true', help='Perform Mendelian Randomization through mixed linear model')
parser_mr.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for Mendelian Randomization analysis')
parser_mr.add_argument('-pvalue', default=1, type=float, metavar='[default:1]', help='The pvalue cutoff of Mendelian randomization analysis result output,default 1')
parser_mr.add_argument('-net', action='store_true', help='Mendelian Randomization network analysis, used to identify functional modules')
parser_mr.add_argument('-o', default='MODAS_mr_out', metavar='[default:MODAS_mr_out]', help='The prefix of the output file, default MODAS_mr_out')
parser_visual = subparsers.add_parser('visual', help='whole genome-wide association analysis and visualization', usage='%(prog)s [options]')
parser_visual.add_argument('-g', metavar='', help='Genotype file in plink bed format for whole genome-wide association analysis')
parser_visual.add_argument('-phe', metavar='', help='Phenotype file for GWAS analysis and visualization')
parser_visual.add_argument('-gwas_model', default='gemma_MLM', type=str, metavar='[default:gemma_MLM]', help='model for GWAS analysis, supports the models provided by gemma, rMVP, and GAPIT software, including gemma_LM, gemma_MLM, rMVP_GLM, rMVP_MLM, rMVP_FarmCPU, GAPIT_GLM, GAPIT_MLM, GAPIT_CMLM, GAPIT_MLMM, GAPIT_FarmCPU, GAPIT_Blink. Default is gemma_MLM.')
parser_visual.add_argument('-qtl', metavar='', help='Phenotype QTL file generated by subcommand regiongwas')
parser_visual.add_argument('-visual', action='store_true', help='Generate web-based GWAS visualization results')
parser_visual.add_argument('-anno', metavar='', help='Gene annotation file for QTL annotation, used to display QTL information in visualized results')
parser_visual.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads genome-wide association analysis and visualization')
parser_visual.add_argument('-o', default='MODAS_visual_out',metavar='[default:MODAS_visual_out]', help='The prefix of the output file,default MODAS_visual_out')
parser_contrast = subparsers.add_parser('contrast', help='construction of stress response index and identification of stress-responsive QTLs', usage='%(prog)s [options]')
parser_contrast.add_argument('-stress_phe', metavar='', help='omics dataset under stress treatment')
parser_contrast.add_argument('-control_phe', metavar='', help='omics dataset under control treatment')
parser_contrast.add_argument('-alpha', default=1, type=int, metavar='[default:1]', help='Contrast strength for contrastive PCA, it represents the trade-off between having the high stress variance and the low control variance, default 1')
parser_contrast.add_argument('-n_comp', default=2, type=int, metavar='[default:2]', help='Number of components for contrastive PCA, default 2')
#parser_contrast.add_argument('-cpca_phe', metavar='', help='stress response index file generated by contrastive PCA, it used to identify stress-responsive molecular QTL')
parser_contrast.add_argument('-beta_test_pvalue', default=None, type=float, metavar='[default:1/number of CPCA PC for beta test]', help='Significance threshold for stress-responsive CPCA PC, used to determine the CPCA PC representing stress response effect of molecular features, default 1/number of CPCA PC for beta test')
parser_contrast.add_argument('-gwas', action='store_true', help='perform GWAS on stress response index to identify stress-responsive molecular QTL')
parser_contrast.add_argument('-g', metavar='', help='Genotype file in plink bed format, used to perform GWAS on stress response index')
parser_contrast.add_argument('-genome_cluster', metavar='', help='Pseudo-genotype file for phenotype pre-screening, generated by subcommand genoidx')
parser_contrast.add_argument('-anova_pvalue', default=0.05, type=float, metavar='[default:0.05]', help='Significance threshold of two-way ANOVA, used to determine stress-responsvie molecular QTL, default 0.05')
parser_contrast.add_argument('-gwas_model', default='MLM', type=str, metavar='[default:MLM]', help='GWAS model for  association analysis, supporting LM, GLM and MLM, default MLM.')
parser_contrast.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for gwas analysis, default 1')
parser_contrast.add_argument('-o', default='MODAS_contrast_out', metavar='[default:MODAS_contrast_out]', help='The prefix of the output file, default MODAS_contrast_out')
parser_coloc = subparsers.add_parser('coloc', help='co-localization analysis of QTLs for multiple traits', usage='%(prog)s [options]')
parser_coloc.add_argument('-qtl', metavar='', help='molecular QTL file for co-localization analysis')
parser_coloc.add_argument('-g', metavar='', help='Genotype file in plink bed format, used to extract genotype of SNPs associated with molecular features')
parser_coloc.add_argument('-gwas_dir', metavar='', help='directory containing GWAS results for molecular features, used to extract SNPs associated with molecular features.It is recommended to use the association analysis results generated by the MODAS subcommand regiongwas')
parser_coloc.add_argument('-metric', default='calinski_harabasz', type=str, metavar='[default:calinski_harabasz]', help='metric for determining clustering patterns of kinship, it can be calinski_harabasz or silhouette, , defualt calinski_harabasz')
parser_coloc.add_argument('-eps', default=0.2, type=float, metavar='[default:0.2]', help='The maximum distance between two traits used to consider a trait to co-localize with the other trait, defualt 0.2')
parser_coloc.add_argument('-pvalue', default=1e-6, type=float, metavar='[defalut:1e-6]', help='Threshold of GWAS results for screening SNPs associated with molecular features, defualt 1e-6')
parser_coloc.add_argument('-p', default=1, type=int, metavar='[default:1]', help='Number of threads for co-localization analysis, default 1')
parser_coloc.add_argument('-o', default='MODAS_coloc_out', metavar='[default:MODAS_coloc_out]', help='The prefix of the output file, default MODAS_coloc_out')


if __name__ == "__main__":
    main()
