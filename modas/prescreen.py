import pandas as pd
import numpy as np
import modas.multiprocess as mp
from sklearn.preprocessing import MinMaxScaler
import os, glob
import logging, re
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.embedded import RRuntimeError
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import subprocess


pandas2ri.activate()
rpy2_logger.setLevel(logging.ERROR)
base = importr('base')
utils = importr('utils')

utils_path = subprocess.check_output('locate modas/utils', shell=True, text=True, encoding='utf-8')
utils_path = '/'.join(re.search('\n(.*site-packages.*)\n', utils_path).group(1).split('/')[:-1])

if not base.require('rMVP')[0]:
    utils.install_packages(np.array(['data.table', 'ggplot2', 'ggsignif', 'Matrix', 'bigmemory', 'RcppProgress', 'BH']), repos='https://cloud.r-project.org', quiet=True)
    utils.install_packages(utils_path + '/rMVP_1.0.6_modify.tar.gz', repos=robjects.rinterface.NULL, type='source', quiet=True)
    utils.install_packages('bigsnpr', dependence=True, repos='https://cloud.r-project.org', quiet=True)
rMVP = importr('rMVP')
bigmemory = importr('bigmemory')


def qtl_pc2bimbam(qtl_pc):
    #qtl_pc.loc[:,:] = np.around(MinMaxScaler(feature_range=(0, 2)).fit_transform(qtl_pc.values),decimals=3)
    g = qtl_pc.T.reset_index()
    g.insert(1,'minor',['A']*g.shape[0])
    g.insert(2,'major',['T']*g.shape[0])
    a = g.iloc[:,0].to_frame()
    a['pos'] = a.iloc[:,0].apply(lambda x: str((int(x.split('_')[2])+int(x.split('_')[3]))/2))
    a['chr'] = a.iloc[:,0].apply(lambda x: x.split('_')[1])
    return a,g


def qtl_pc_gwas_parallel(omics_phe, bimbam_dir, threads, geno, geno_prefix, gwas_model):
    qtl_pc_gwas_args = list()
    if gwas_model == 'MLM' or gwas_model == 'GLM':
        fam = pd.read_csv(geno + '.fam', sep=r'\s+', header=None)
        fam[5] = 1
        fam.to_csv(geno_prefix + '.link.fam', sep='\t', na_rep='NA', header=None, index=False)
        #omics_phe = omics_phe.reindex(fam[0].values)
        if os.path.exists(geno_prefix + '.link.bed'):
            os.remove(geno_prefix + '.link.bed')
        if os.path.exists(geno_prefix + '.link.bim'):
            os.remove(geno_prefix + '.link.bim')
        os.symlink(geno + '.bed', geno_prefix + '.link.bed')
        os.symlink(geno + '.bim', geno_prefix + '.link.bim')
        if gwas_model == 'MLM':
            related_matrix_cmd = utils_path + '/gemma -bfile {0}.link -gk 1 -o {1}'.format(geno_prefix,geno_prefix)
            s = mp.run(related_matrix_cmd)
            if s!=0:
                return None
    if gwas_model == 'MLM':
        gemma_cmd = utils_path + '/gemma -g {0} -a {1} -p {2} -k ./output/{3}.cXX.txt -lmm -n 1 -o {4}'
    elif gwas_model == 'LM':
        gemma_cmd = utils_path + '/gemma -g {0} -a {1} -p {2} -lm  -o {3}'
    else:
        g = pd.read_csv(bimbam_dir.strip('/')+'/'+geno_prefix+'_qtl_pc.geno.txt',header=None)
        a = pd.read_csv(bimbam_dir.strip('/')+'/'+geno_prefix+'_qtl_pc.anno.txt',header=None)
        g.iloc[:,3:].to_csv(geno_prefix+'.numeric.txt',index=False, header=None, sep='\t')
        a.columns = ['SNP', 'Pos', 'Chr']
        a = a[['SNP', 'Chr', 'Pos']]
        a.to_csv(geno_prefix+'.map.txt', index=False, sep='\t')

    for m in omics_phe.columns:
        phe = omics_phe[m].to_frame()
        m = m.replace('m/z', 'm.z')
        phe.to_csv(bimbam_dir.strip('/') + '/' + m + '_phe.txt', index=False, header=None, na_rep='NA')
        if gwas_model == 'MLM':
            qtl_pc_gwas_args.append((gemma_cmd.format(bimbam_dir.strip('/') + '/'+geno_prefix+'_qtl_pc.geno.txt', bimbam_dir.strip('/') + '/'+geno_prefix+'_qtl_pc.anno.txt', bimbam_dir.strip('/') + '/' + m + '_phe.txt',geno_prefix, m + '_prescreen'),))
        elif gwas_model == 'LM':
            qtl_pc_gwas_args.append((gemma_cmd.format(bimbam_dir.strip('/')+'/'+geno_prefix+'_qtl_pc.geno.txt',bimbam_dir.strip('/')+'/'+geno_prefix+'_qtl_pc.anno.txt',bimbam_dir.strip('/')+'/'+m+'_phe.txt',m+'_prescreen'),))
        #else:
        #    qtl_pc_gwas_args.append((phe.reset_index(), geno_prefix+'.link', geno_prefix+'.numeric.txt', geno_prefix+'.map.txt', 1, './output'))
    if gwas_model == 'LM' or gwas_model == 'MLM':
        s = mp.parallel(mp.run, qtl_pc_gwas_args, threads)
    else:
        if not os.path.exists('./output'):
            os.mkdir('./output')
        #s = mp.parallel(glm_gwas, (qtl_pc_gwas_args[0],), 1)
        #s = mp.parallel(glm_gwas, qtl_pc_gwas_args[1:], threads)
        omics_phe.columns = [i.replace('m/z', 'm.z') for i in omics_phe.columns]
        s = glm_gwas(omics_phe, geno_prefix+'.link', geno_prefix+'.numeric.txt', geno_prefix+'.map.txt', 1, './output')
    if gwas_model == 'MLM' or gwas_model == 'GLM':
        os.remove(geno_prefix+'.link.bed')
        os.remove(geno_prefix+'.link.bim')
        os.remove(geno_prefix+'.link.fam')
    return s


def glm_gwas(omics_phe, pc_geno_prefix, genofile, mapfile, threads, out_path):
    try:
        geno_prefix = '.'.join(genofile.split('/')[-1].split('.')[:-2])
        base.sink('/dev/null')
        robjects.r('''
                gwas <- function(omics_phe, pc_geno_prefix, geno_prefix, genofile, mapfile, threads, out_path){
                    library(rMVP)
                    if(!file.exists(paste(pc_geno_prefix,'.pc.desc',sep=''))){
                         MVP.Data(fileBed=pc_geno_prefix, fileKin=F, filePC=F, out=pc_geno_prefix, verbose=F)
                         MVP.Data.PC(T,mvp_prefix=pc_geno_prefix, pcs.keep=5, verbose=F)
                    }
                    MVP.Data(fileNum=genofile, fileMap=mapfile, fileKin=F, filePC=F, sep_num='\t', type.geno='double', out=geno_prefix)
                    geno = attach.big.matrix(paste(geno_prefix, '.geno.desc',sep=''))
                    map_file = read.table(paste(geno_prefix, '.geno.map',sep=''),sep='\t',header=T)
                    Covariates_PC = bigmemory::as.matrix(attach.big.matrix(paste(pc_geno_prefix,'.pc.desc',sep='')))
                    phe_name = names(omics_phe)
                    for(i in 2:ncol(omics_phe)){
                        mvp = MVP(phe=omics_phe[,c(1,i)], geno=geno, map=map_file, CV.GLM=Covariates_PC, priority='speed', nPC.GLM=5,
                        ncpus=threads, maxLoop=10, threshold=0.05, method=c('GLM'), file.output=F, verbose=F)
                        res = cbind(mvp$map, mvp$glm.results)
                        names(res) <- c('rs', 'chr', 'ps', 'effect', 'se', 'p_wald')
                        print(head(res))
                        write.table(res,file=paste(out_path,'/',as.character(phe_name[i]),'_prescreen.assoc.txt',sep=''),sep='\t', quote=F, row.names=F)
                    }
                }
                ''')
        gwas = robjects.r('gwas')
        gwas(omics_phe, pc_geno_prefix, geno_prefix, genofile, mapfile, threads, out_path)
        base.sink()
    except RRuntimeError:
        return 0
    except ValueError:
        return 0
    else:
        return 1

# def glm_gwas(omics_phe, pc_geno_prefix, genofile, mapfile, threads, out_path):
#     try:
#         geno_prefix = '.'.join(genofile.split('/')[-1].split('.')[:-2])
#         base.sink('/dev/null')
#         if not os.path.exists(pc_geno_prefix + '.pc.desc'):
#             rMVP.MVP_Data(fileBed=pc_geno_prefix, fileKin=False, filePC=False, out=pc_geno_prefix,
#                           verbose=False)
#             rMVP.MVP_Data_PC(True, mvp_prefix=pc_geno_prefix, pcs_keep=5, verbose=False)
#         rMVP.MVP_Data(fileNum=genofile, fileMap=mapfile, fileKin=False, filePC=False, sep_num = '\t', type_geno='double',out=geno_prefix, verbose=False)
#         geno = bigmemory.attach_big_matrix(geno_prefix +'.geno.desc')
#         map_file = pd.read_csv(geno_prefix +'.geno.map', sep='\t')
#         Covariates_PC = bigmemory.as_matrix(bigmemory.attach_big_matrix(pc_geno_prefix + '.pc.desc'))
#         # base.setwd('./output')
#         mvp = rMVP.MVP(phe=omics_phe, geno=geno, map=map_file, CV_GLM=Covariates_PC, priority="speed", nPC_GLM=5,
#                 ncpus=threads, maxLoop=10, threshold=0.05, method=['GLM'], file_output=False, verbose=False)
#         gwas_res = pd.DataFrame(mvp.rx2('glm.results'), columns=['effect', 'se', 'p_wald'])
#         pos = pd.DataFrame(mvp.rx2('map'))
#         pos.columns = ['rs','chr', 'ps']
#         pos.index = gwas_res.index
#         res = pd.concat([pos, gwas_res], axis=1)
#         res.to_csv(out_path.rstrip('/') + '/' + str(omics_phe.columns[1]) + '_prescreen.assoc.txt', index=False,sep='\t')
#         base.sink()
#     except RRuntimeError:
#         return 0
#     except ValueError:
#         return 0
#     else:
#         return 1


def prescreen(omics_phe,suggest_pvalue):
    phe_sig_qtl = list()
    sig_phe_names = list()
    for fn in glob.glob('output/*_prescreen.assoc.txt'):
        gwas = pd.read_csv(fn, sep='\t')
        if gwas['p_wald'].min() > suggest_pvalue:
            continue
        phe_name = fn.split('/')[-1].replace('_prescreen.assoc.txt','')
        pos = list()
        for rs in gwas.loc[gwas.p_wald <= suggest_pvalue, 'rs'].values:
            chrom,start,end = rs.split('_')[1:4]
            start,end = int(start),int(end)
            if not pos:
                pos.append([chrom,start,end,phe_name])
            else:
                if pos[-1][0] != chrom:
                    pos.append([chrom,start,end,phe_name])
                else:
                    if start < pos[-1][2]:
                        start = start if start < pos[-1][1] else pos[-1][1]
                        end = end if end > pos[-1][2] else pos[-1][2]
                        pos[-1] = [chrom,start,end,phe_name]
                    else:
                        pos.append([chrom,start,end,phe_name])
        phe_sig_qtl.extend(pos)
        if not gwas.loc[gwas.p_wald <= suggest_pvalue,:].empty:
            sig_phe_names.append(phe_name.replace('m.z','m/z'))
    sig_omics_phe = omics_phe.loc[:, sig_phe_names]
    phe_sig_qtl = pd.DataFrame(phe_sig_qtl,columns=['chr','start','end','phe_name'])
    return sig_omics_phe, phe_sig_qtl



