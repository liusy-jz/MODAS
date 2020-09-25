#!/usr/bin/env python
import sys, time
import modas.genoidx as gi
import modas.phenorm as pn
import modas.prescreen as ps
import modas.localgwas as lg
import pandas as pd
import numpy as np
import datetime
import argparse
import traceback
import warnings
import shutil
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

def genoidx(args,log):
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
        log.log('snp clumping is done.')
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
        genome_cluster = gi.genome_cluster(g_clump,args.w,args.s,args.p)
        if genome_cluster.empty:
            log.log('genome cluster analysis failed,Please check genotype file.')
            sys.exit()
        genome_cluster.to_csv(args.o+'.genome_cluster.csv')
    log.log('genoidx is finishd!')

def phenorm(args,log):
    log.log('begin phenotype filter/normalize analysis...')
    if not os.path.exists(args.phe):
        raise ArgsError('the phenotype file {} is not exists.'.format(args.phe))
    phe = pd.read_csv(args.phe, index_col=0)
    if phe.empty:
        log.log('the phenotype file {} is empty and software is terminated.'.format(args.phe))
        sys.exit()
    if args.v is None and args.r is None and not args.log2 and not args.log10 and not args.ln and not args.norm:
        log.log('No phenotype filter process or normalize process choosed, please select one of the argument in -a, -r, -log2, -log10, -ln, -norm.')
        sys.exit()
    log.log('Read phenotype file finished,There are {} phenotypes in phenotype file.'.format(phe.shape[1]))
    phe = phe.loc[:,(phe.var()-0) >= 1e-6]
    if args.v is not None:
        if args.v <= 0:
            raise ArgsError('metabolite abundance cutoff should larger than zero.')
        log.log('filter phenotype values less than {}'.format(args.v))
        phe = pn.abundance_filter(phe, args.v)
        log.log('after filtered, There are {} phenotypes left.'.format(phe.shape[1]))
    if args.r is not None:
        if args.r < 0 or args.r > 1:
            raise ArgsError('the interval of missing ration is between 0 and 1.')
        log.log('filter phenotype missing ratio larger than {}'.format(args.r))
        phe = pn.missing_filter(phe, args.r)
        log.log('after filtered, There are {} phenotypes left.'.format(phe.shape[1]))
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
    phe.to_csv(args.o+'.normalized_phe.csv',na_rep=0)

def prescreen(args,log):
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
    if not os.path.exists(args.g+'.bed'):
        raise ArgsError('the plink bed format genotype file {} is not exists.'.format(args.g))
    genome_cluster = pd.read_csv(args.genome_cluster, index_col=0)
    if genome_cluster.empty:
        log.log('genome cluster file {} is empty and software is terminated.'.format(args.genome_cluster))
        sys.exit()
    if os.path.exists('tmp_metabolite_bimbam'):
        shutil.rmtree('tmp_metabolite_bimbam')
    os.mkdir('tmp_metabolite_bimbam')
    if os.path.exists('output'):
        os.rename('output','output'+datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))
    fam = pd.read_csv(args.g+'.fam', sep=r'\s+', header=None)
    genome_cluster = genome_cluster.reindex(fam[0].values)
    phe = phe.reindex(fam[0].values)
    log.log('convert genome cluster file to bimbam genotype file.')
    a,g = ps.qtl_pc2bimbam(genome_cluster)
    a.to_csv('tmp_metabolite_bimbam/'+args.g.split('/')[-1]+'_qtl_pc.anno.txt',index=False,header=None)
    g.to_csv('tmp_metabolite_bimbam/'+args.g.split('/')[-1]+'_qtl_pc.geno.txt',index=False,header=None)
    if args.lm_suggest_pvalue is None:
        args.lm_suggest_pvalue = 1.0 / genome_cluster.shape[1]
    log.log('run lm model for genome cluster pseudo-SNP filter.')
    s = ps.qtl_pc_lm_gwas_parallel(phe,'tmp_metabolite_bimbam',args.p,args.g)
    ps.generate_metabolite_qtl_pc_bimbam(phe,a,g,args.lm_suggest_pvalue,args.p)
    log.log('run lmm model for Significant phenotype filter.')
    s = ps.qtl_pc_lmm_gwas_parallel(phe,'tmp_metabolite_bimbam',args.p,args.g)
    sig_omics_phe, phe_sig_qtl = ps.prescreen(phe,args.lmm_suggest_pvalue)
    sig_omics_phe.to_csv(args.o+'.sig_omics_phe.csv')
    phe_sig_qtl.to_csv(args.o+'.phe_sig_qtl.csv',index=False)
    log.log('prescreened phenotype saved successfully')

def localgwas(args,log):
    log.log('begin localgwas analysis...')
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
    lg.generate_qtl_batch(phe,phe_sig_qtl,args.g,args.p,'tmp_qtl_bed','tmp_rs_dir')
    log.log('begin local gwas analysis of phenotype significant gwas region...')
    s = lg.local_gwas_parallel('tmp_qtl_bed', args.p,args.g)
    log.log('begin generate QTL of phenotypes...')
    s = lg.generate_clump_input('output',args.p)
    s = lg.plink_clump('tmp_qtl_bed', args.p1, args.p2, args.p)
    qtl_res = lg.generate_qtl('clump_result')
    qtl_res.to_csv(args.o+'.local_gwas_qtl_res.csv',index=False)
    log.log('generate QTL done.')
    log.log('local gwas analysis is finished.')


def MODAS(args):
    if args.command == 'genoidx':
        log = Logger(args.o+'.genoidx.log.txt')
        genoidx(args, log)
    if args.command == 'phenorm':
        log = Logger(args.o+'.phenorm.log.txt')
        phenorm(args, log)
    if args.command == 'prescreen':
        log = Logger(args.o+'.prescreen.log.txt')
        prescreen(args, log)
    if args.command == 'localgwas':
        log = Logger(args.o+'.localgwas.log.txt')
        localgwas(args, log)


USAGE = ' MODAS: Multi-omics data association study '

parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)
parser_genoidx = subparsers.add_parser('genoidx', help='generate a pseudo-genotype file', usage='%(prog)s [options]')
parser_genoidx.add_argument('-g', metavar='', help='input genotype file(hapmap/plink bed format), plink genotype only need bed format prefix')
parser_genoidx.add_argument('-w', default=1000000, type=float, metavar='[default:1000000]', help='window size for genome cluster analysis')
parser_genoidx.add_argument('-s', default=500000, type=float, metavar='[default:500000]', help='step size for genome cluster analysis')
parser_genoidx.add_argument('-genome_cluster', action='store_true', help='genome cluster analysis for input genotype file')
parser_genoidx.add_argument('-convert', action='store_true', help='convert hapmap genotype format to plink bed format genotype')
parser_genoidx.add_argument('-clump', action='store_true', help='clumping analysis for genotype file')
parser_genoidx.add_argument('-p', default=1, type=int, metavar='[default:1]', help='number of threads for genome cluster analysis,max threads is genome chromosome number.')
parser_genoidx.add_argument('-o', default='MODAS_genoidx_output', metavar='', help='output file prefix')
parser_phenorm = subparsers.add_parser('phenorm', help='phenotype normalization and transformation', usage='%(prog)s [options]')
parser_phenorm.add_argument('-phe', metavar='', help='input phenotype file')
parser_phenorm.add_argument('-r', default=None, type=float, metavar='', help='phenotype missing ratio filter cutoff')
parser_phenorm.add_argument('-v', default=None, type=float, metavar='', help='phenotype value filter cutoff')
parser_phenorm.add_argument('-log2', action='store_true', help='log2 transform phenotype')
parser_phenorm.add_argument('-log10', action='store_true', help='log10 trasform phenotype')
parser_phenorm.add_argument('-ln', action='store_true', help='ln transform phenotype')
parser_phenorm.add_argument('-norm', action='store_true', help='normalize phenotype by boxcox method')
parser_phenorm.add_argument('-o', default='MODAS_phenorm_output', metavar='', help='output file prefix')
parser_prescreen = subparsers.add_parser('prescreen', help='prescreening of associated QTLs with lm/lmm', usage='%(prog)s [options]')
parser_prescreen.add_argument('-g', metavar='', help='input plink bed genotype file')
parser_prescreen.add_argument('-genome_cluster', metavar='', help='genome cluster file generated by genoidx sub-command')
parser_prescreen.add_argument('-phe', metavar='', help='input phenotype file')
parser_prescreen.add_argument('-lm_suggest_pvalue', default=None, type=float, metavar='[default:1/ genome cluster file converted pseudo-SNP number]', help='lm suggest pvalue for pseudo-SNP filter' )
parser_prescreen.add_argument('-lmm_suggest_pvalue', default=1e-6, type=float, metavar='[default:1e-6]', help='lmm suggest pvalue for phenotype prescreen')
parser_prescreen.add_argument('-p', default=1, type=int, metavar='[default:1]', help='number of threads for analysis in prescreen sub command')
parser_prescreen.add_argument('-o', default='MODAS_prescreen_output', metavar='', help='output file prefix')
parser_localgwas = subparsers.add_parser('localgwas', help='perform local gwas to identify QTLs', usage='%(prog)s [options]')
parser_localgwas.add_argument('-g', metavar='', help='input plink bed genotype file')
parser_localgwas.add_argument('-phe', metavar='', help='input phenotype file')
parser_localgwas.add_argument('-phe_sig_qtl', metavar='', help='significant phenotype QTL region file')
parser_localgwas.add_argument('-p1', metavar='', help='significance threshold for index SNPs')
parser_localgwas.add_argument('-p2', metavar='', help='secondary significance threshold for clumped SNPs')
parser_localgwas.add_argument('-p', default=1, type=int, metavar='[default:1]', help='number of threads for local gwas analysis')
parser_localgwas.add_argument('-o', default='MODAS_localgwas_out', metavar='', help='output file prefix')

if len(sys.argv) <= 2:
    sys.argv.append('-h')
args = parser.parse_args()
MODAS(args)
