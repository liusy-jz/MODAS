import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.embedded import RRuntimeError
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import subprocess
import logging
import glob, os
import shutil
import re

pandas2ri.activate()
rpy2_logger.setLevel(logging.ERROR)
rMVP = importr('rMVP')
base = importr('base')
data_table = importr('data.table')
bigmemory = importr('bigmemory')


utils_path = subprocess.check_output('locate modas/utils', shell=True, text=True, encoding='utf-8')
utils_path = '/'.join(re.search('\n(.*site-packages.*)\n', utils_path).group(1).split('/')[:-1])


def gemma_cmd(model, geno_prefix, kin_prefix, n, out_prefix):
    if model == 'LM':
        return utils_path + '/gemma -bfile {0} -lm  -o {1}'.format(geno_prefix, out_prefix)
    if model == 'MLM':
        return utils_path + '/gemma -bfile {0} -k ./output/{1}.cXX.txt -lmm -n {2} -o {3}'.format(geno_prefix, kin_prefix, n, out_prefix)


def rmvp(model, cv_geno_prefix, geno_prefix, omics_phe, threads, out_path):
    try:
        base.sink('/dev/null')
        if model == 'GLM' or model == 'FarmCPU':
            if not os.path.exists(cv_geno_prefix + '.pc.desc'):
                rMVP.MVP_Data(fileBed=cv_geno_prefix, fileKin=False, filePC=False, out=cv_geno_prefix,
                              verbose=False)
                rMVP.MVP_Data_PC(True, mvp_prefix=cv_geno_prefix, pcs_keep=10, verbose=False)
        if model == 'MLM':
            if not os.path.exists(cv_geno_prefix + '.kin.desc'):
                rMVP.MVP_Data(fileBed=cv_geno_prefix, fileKin=False, filePC=False, out=cv_geno_prefix,
                              verbose=False)
                rMVP.MVP_Data_Kin(True, mvp_prefix=cv_geno_prefix, verbose=False)
        if not os.path.exists(geno_prefix + '.geno.desc'):
            rMVP.MVP_Data(fileBed=geno_prefix, fileKin=False, filePC=False, out=geno_prefix, verbose=False)
        geno = bigmemory.attach_big_matrix(geno_prefix +'.geno.desc')
        map_file = pd.read_csv(geno_prefix +'.geno.map', sep='\t')
        if model == 'GLM' or model == 'FarmCPU':
            Covariates_PC = bigmemory.as_matrix(bigmemory.attach_big_matrix(cv_geno_prefix + '.pc.desc'))
        if model == 'MLM':
            Kinship = bigmemory.attach_big_matrix(cv_geno_prefix + '.kin.desc')
        if model == 'GLM':
            # robjects.r('''
            #     gwas <- function(omics_phe, geno, map_file, Covariates_PC, threads){
            #         library(rMVP)
            #         mvp <- MVP(phe=omics_phe, geno=geno, map=map_file, CV.GLM=Covariates_PC, priority='speed', nPC.GLM=5,
            #         ncpus=threads, maxLoop=10, threshold=0.05, method=c('GLM'), file.output=F, verbose=F)
            #         res <- cbind(mvp$map, mvp$glm.results)
            #         return(res)
            #     }
            # ''')
            # mvp = robjects.r('gwas')
            # res = mvp(omics_phe, geno, map_file, Covariates_PC, threads)
            mvp = rMVP.MVP(phe=omics_phe, geno=geno, map=map_file, CV_GLM=Covariates_PC, priority='speed', nPC_GLM=5,
                    ncpus=threads, maxLoop=10, threshold=0.05, method=['GLM'], file_output=False,
                    verbose=False)
            gwas_res = pd.DataFrame(mvp.rx2('glm.results'), columns=['Effect', 'SE', str(omics_phe.columns[1]) + '.GLM'])
            pos = pd.DataFrame(mvp.rx2('map'))
            pos.index = gwas_res.index
            res = pd.concat([pos, gwas_res], axis=1)
        if model == 'FarmCPU':
            # robjects.r('''
            #     gwas <- function(omics_phe, geno, map_file, Covariates_PC, threads){
            #         library(rMVP)
            #         mvp <- MVP(phe=omics_phe, geno=geno, map=map_file, CV.GLM=Covariates_PC, priority='speed', nPC.GLM=5,
            #         ncpus=threads, maxLoop=10, threshold=0.05, method=c('FarmCPU'), method.bin='static', file.output=F, verbose=F)
            #         res <- cbind(mvp$map, mvp$farmcpu.results)
            #         return(res)
            #     }
            # ''')
            # mvp = robjects.r('gwas')
            # res = mvp(omics_phe, geno, map_file, Covariates_PC, threads)
            mvp = rMVP.MVP(phe=omics_phe, geno=geno, map=map_file, CV_FarmCPU=Covariates_PC, priority='speed', nPC_FarmCPU=3,
                     ncpus=threads, maxLoop=10, threshold=0.05, method=['FarmCPU'], file_output=False, method_bin='static',
                     verbose=True)
            gwas_res = pd.DataFrame(mvp.rx2('farmcpu.results'), columns=['Effect', 'SE', str(omics_phe.columns[1]) + '.FarmCPU'])
            pos = pd.DataFrame(mvp.rx2('map'))
            pos.index = gwas_res.index
            res = pd.concat([pos, gwas_res], axis=1)
        if model == 'MLM':
            # robjects.r('''
            #     gwas <- function(omics_phe, geno, map_file, Kinship, threads){
            #         library(rMVP)
            #         mvp <- MVP(phe=omics_phe, geno=geno, map=map_file, K=Kinship, priority='speed', nPC.GLM=5,
            #         vc.method='BRENT', ncpus=threads, maxLoop=10, threshold=0.05, method=c('MLM'), file.output=F, verbose=F)
            #         res <- cbind(mvp$map, mvp$mlm.results)
            #         return(res)
            #     }
            # ''')
            # mvp = robjects.r('gwas')
            # res = mvp(omics_phe, geno, map_file, Kinship, threads)
            mvp = rMVP.MVP(phe=omics_phe, geno=geno, map=map_file, K=Kinship, priority='speed', vc_method='BRENT',
                     ncpus=threads, maxLoop=10, threshold=0.05, method=['MLM'], file_output=False,
                     verbose=False)
            gwas_res = pd.DataFrame(mvp.rx2('mlm.results'), columns=['Effect', 'SE', str(omics_phe.columns[1])+'.MLM'])
            pos = pd.DataFrame(mvp.rx2('map'))
            pos.index = gwas_res.index
            res = pd.concat([pos, gwas_res], axis=1)
        res.to_csv(out_path.rstrip('/') + '/' + str(omics_phe.columns[1])+'.' + model + '.csv', index=False)
        base.sink()
    except RRuntimeError:
        return 1
    except ValueError:
        return 1
    else:
        return 0


def gapit(model, geno, omics_phe, gapit_path):
    try:
        base.sink('/dev/null')
        robjects.r('source("'+gapit_path.rstrip('/')+'/GAPIT.library.R")')
        robjects.r('source("'+gapit_path.rstrip('/')+'/gapit_functions.txt")')
        robjects.r('''gapit <- function(geno,omics_phe,model){
            library(bigsnpr)
            g <- snp_readBed(paste(geno,'.bed',sep=''), backingfile=tempfile())
            g <- snp_attach(g)
            GD <- cbind(g$fam$family.ID,as.data.frame(snp_fastImputeSimple(g$genotypes, method='mode')[]))
            names(GD) <- c('Taxa',g$map$marker.ID)
            GM <- g$map[c('marker.ID','chromosome','physical.pos')]
            names(GM) <- c('Name','Chromosome','Position')
            GAPIT(Y=omics_phe, GD=GD, GM=GM, model=model, Major.allele.zero = T, SNP.MAF=0.05)
        }''')
        GAPIT = robjects.r('gapit')
        GAPIT(geno, omics_phe, model)
        base.sink()
    except RRuntimeError:
        return 1
    except ValueError:
        return 1
    else:
        return 0
