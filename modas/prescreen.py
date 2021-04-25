import pandas as pd
import numpy as np
import modas.multiprocess as mp
from sklearn.preprocessing import MinMaxScaler
import os,glob

def qtl_pc2bimbam(qtl_pc):
    #qtl_pc.loc[:,:] = np.around(MinMaxScaler(feature_range=(0, 2)).fit_transform(qtl_pc.values),decimals=3)
    g = qtl_pc.T.reset_index()
    g.insert(1,'minor',['A']*g.shape[0])
    g.insert(2,'major',['T']*g.shape[0])
    a = g.iloc[:,0].to_frame()
    a['pos'] = a.iloc[:,0].apply(lambda x: str((int(x.split('_')[2])+int(x.split('_')[3]))/2))
    a['chr'] = a.iloc[:,0].apply(lambda x: x.split('_')[1])
    return a,g


def qtl_pc_genotype_subset(g,a,rs,output_dir,phe):
    sub_g = g.loc[g.iloc[:,0].isin(rs),:]
    sub_a = a.loc[g.iloc[:,0].isin(rs),:]
    sub_g.to_csv(output_dir+'tmp_'+phe.replace('m/z','m.z')+'.geno.txt',index=False,header=None)
    sub_a.to_csv(output_dir+'tmp_'+phe.replace('m/z','m.z')+'.anno.txt',index=False,header=None)


def generate_omics_qtl_pc_bimbam(omics_phe,a,g,lm_suggest_pvalue,threads):
    prefix = 'tmp_omics_phe_bimbam/'
    #subset_genotype_args = list()
    g.index = g['index']
    a.index = g['index']
    for phe in omics_phe.columns:
        gwas = pd.read_csv('output/'+phe.replace('m/z','m.z')+'_bimbam_lm.assoc.txt',sep='\t')
        rs = gwas.loc[gwas.p_wald <= lm_suggest_pvalue,'rs']
        #index = g.iloc[:,0].isin(rs)
        #sub_g = g.loc[index,:]
        #sub_a = a.loc[index,:]
        sub_g = g.reindex(rs)
        sub_a = a.reindex(rs)
        sub_g.to_csv(prefix+'tmp_'+phe.replace('m/z','m.z')+'.geno.txt',index=False,header=None)
        sub_a.to_csv(prefix+'tmp_'+phe.replace('m/z','m.z')+'.anno.txt',index=False,header=None)
#        subset_genotype_args.append((g,a,rs,prefix,phe))
#    s = mp.parallel(qtl_pc_genotype_subset,subset_genotype_args, threads)


def qtl_pc_lm_gwas_parallel(omics_phe,bimbam_dir,threads,geno):
    qtl_pc_lm_args = list()
    geno_prefix = geno.split('/')[-1]
    gemma_cmd = 'gemma.linux -g {0} -a {1} -p {2} -lm  -o {3}'
    for m in omics_phe.columns:
        phe = omics_phe[m].to_frame()
        m = m.replace('m/z','m.z')
        phe.to_csv(bimbam_dir.strip('/')+'/'+m+'_phe.txt',index=False,header=None,na_rep='NA')
        qtl_pc_lm_args.append((gemma_cmd.format(bimbam_dir.strip('/')+'/'+geno_prefix+'_qtl_pc.geno.txt',bimbam_dir.strip('/')+'/'+geno_prefix+'_qtl_pc.anno.txt',bimbam_dir.strip('/')+'/'+m+'_phe.txt',m+'_bimbam_lm'),))
    s = mp.parallel(mp.run, qtl_pc_lm_args, threads)
    return s


def qtl_pc_lmm_gwas_parallel(omics_phe,bimbam_dir,threads,geno,sample_id):
    qtl_pc_lmm_args = list()
    #g = read_plink1_bin(geno+'.bed', geno+'.bim', geno+'.fam', verbose=False)
    #g = g.sel(sample=sample_id)
    geno_prefix = geno.split('/')[-1]
    #if os.path.exists(geno_prefix+'.link.bed'):
    #    os.remove(geno_prefix+'.link.bed')
    #if os.path.exists(geno_prefix+'.link.bim'):
    #    os.remove(geno_prefix+'.link.bim')
    #write_plink1_bin(g,geno_prefix+'.link.bed', geno_prefix+'.link.bim,', geno_prefix+'.link.fam',verbose=False)
    fam = pd.read_csv(geno+'.fam', sep=r'\s+', header=None)
    fam[5] = 1
    fam.to_csv(geno_prefix+'.link.fam', sep='\t', na_rep='NA', header=None, index=False)
    omics_phe = omics_phe.reindex(fam[0].values)
    omics_phe.to_csv('bimbam_phe.txt',sep='\t',index=False,header=None,na_rep='NA')
    if os.path.exists(geno_prefix+'.link.bed'):
        os.remove(geno_prefix+'.link.bed')
    if os.path.exists(geno_prefix+'.link.bim'):
        os.remove(geno_prefix+'.link.bim')
    os.symlink(geno+'.bed', geno_prefix+'.link.bed')
    os.symlink(geno+'.bim', geno_prefix+'.link.bim')
    related_matrix_cmd = 'gemma.linux -bfile {0}.link -gk 1 -o {1}'.format(geno_prefix,geno_prefix)
    s = mp.run(related_matrix_cmd)
    if s!=0:
        return None
    gemma_cmd = 'gemma.linux -g {0} -a {1} -p bimbam_phe.txt -k ./output/{2}.cXX.txt -lmm -n {3} -o {4}'
    for _,m in enumerate(omics_phe.columns):
        m = m.replace('m/z','m.z')
        qtl_pc_lmm_args.append((gemma_cmd.format(bimbam_dir.strip('/')+'/tmp_'+m+'.geno.txt',bimbam_dir.strip('/')+'/tmp_'+m+'.anno.txt',geno_prefix,_+1,m+'_bimbam'),))
    s = mp.parallel(mp.run, qtl_pc_lmm_args, threads)
    os.remove(geno_prefix+'.link.bed')
    os.remove(geno_prefix+'.link.bim')
    os.remove(geno_prefix+'.link.fam')
    return s

def prescreen(omics_phe,lmm_suggest_pvalue):
    phe_sig_qtl = list()
    sig_phe_names = list()
    for fn in glob.glob('output/*_bimbam.assoc.txt'):
        gwas = pd.read_csv(fn,sep='\t')
        if gwas['p_wald'].min() > lmm_suggest_pvalue:
            continue
        phe_name = fn.split('/')[-1].replace('_bimbam.assoc.txt','')
        pos = list()
        for rs in gwas.loc[gwas.p_wald <= lmm_suggest_pvalue, 'rs'].values:
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
        if not gwas.loc[gwas.p_wald <= lmm_suggest_pvalue,:].empty:
            sig_phe_names.append(phe_name.replace('m.z','m/z'))
    sig_omics_phe = omics_phe.loc[:, sig_phe_names]
    phe_sig_qtl = pd.DataFrame(phe_sig_qtl,columns=['chr','start','end','phe_name'])
    return sig_omics_phe, phe_sig_qtl



