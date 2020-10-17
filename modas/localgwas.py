import pandas as pd
import numpy as np
import modas.multiprocess as mp
from sklearn.decomposition import PCA
import glob, os
import shutil
import re


def local_gwas_parallel(bed_dir, threads,geno):
    local_gwas_args = list()
    geno_prefix = geno.split('/')[-1]
    fam = pd.read_csv(geno+'.fam', sep=r'\s+', header=None)
    fam[5] = 1
    fam.to_csv(geno_prefix+'.link.fam', sep='\t', na_rep='NA', header=None, index=False)
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
    gemma_cmd = 'gemma.linux -bfile {0} -k ./output/{1}.cXX.txt -lmm -n 1 -o {2}'
    for i in glob.glob(bed_dir+'/*.bed'):
        i = i.replace('.bed','')
        i = i.replace('m/z','m.z')
        prefix = i.split('/')[-1]
        local_gwas_args.append((gemma_cmd.format(i,geno_prefix,prefix+'_plink'),))
    s = mp.parallel(mp.run, local_gwas_args, threads)
    os.remove(geno_prefix+'.link.bed')
    os.remove(geno_prefix+'.link.bim')
    os.remove(geno_prefix+'.link.fam')
    return s

def generate_qtl_batch(omics_phe,phe_sig_qtl,geno_name,threads,bed_dir,rs_dir):
    plink_extract = 'plink -bfile {} --extract {} --make-bed -out {}'
    bim = pd.read_csv(geno_name+'.bim',sep='\t',header=None)
    qtl_batch = list()
    rs = dict()
    for index,row in phe_sig_qtl.iterrows():
        rs.setdefault(row['phe_name'],[]).extend(bim.loc[(bim[0]==row['chr']) & (bim[3]>=row['start']) & (bim[3]<=row['end']),1].values.tolist())
    for phe_name in rs:
        out_name = bed_dir.strip('/') + '/' + '_'.join(['tmp',phe_name])
        rs_name = rs_dir.strip('/') + '/' + '_'.join(['tmp',phe_name,'rs.txt'])
        pd.Series(rs[phe_name]).to_frame().to_csv(rs_name,index=False,header=False)
        qtl_batch.append((plink_extract.format(geno_name,rs_name,out_name),))
    mp.parallel(mp.run,qtl_batch,threads)
    for fn in glob.glob(bed_dir.strip('/')+'/*fam'):
        fam = pd.read_csv(fn,sep=' ',header=None)
        phe_name = '_'.join(fn.split('/')[-1].split('_')[1:]).replace('m.z','m/z').replace('.fam','')
        fam.loc[:,5] = omics_phe.loc[:,phe_name].reindex(fam.loc[:,0]).values
        fam.to_csv(fn,index=False,header=None,sep=' ',na_rep='NA')


def generate_clump_input(dir,num_threads):
    if os.path.exists('./clump_input'):
        shutil.rmtree('./clump_input')
    os.mkdir('./clump_input')
    cmd = '''awk '{if(NR==1)print "SNP\\tP"; else print $2"\\t"$12}' '''
    cmds = list()
    fns = list()
    for fn in glob.glob(dir.strip('/')+'/*_plink.assoc.txt'):
        filename = fn.split('/')[-1]
        cmds.append((cmd+'{0} > ./clump_input/{1}'.format(fn, filename.replace('_plink.assoc.txt', '.assoc')),))
        fns.append(filename)
    s = mp.parallel(mp.run, cmds, num_threads)
    if sum(s) != 0:
        print(','.join(list(np.array(fns)[s]))+' do not  successfully generated clump input file.')
    return s


def plink_clump(geno_path, p1, p2, num_threads):
    if os.path.exists('./clump_result'):
        shutil.rmtree('./clump_result')
    os.mkdir('./clump_result')
    cmd = 'plink --bfile {0} --clump {1}  --clump-p1 {2} --clump-p2 {3} --clump-kb {4} --clump-r2 0.2 --out {5}'
    cmds = list()
    ms = list()
    for fn in glob.glob('./clump_input/*'):
        phe_name = fn.split('/')[-1].replace('.assoc','')
        cmds.append((cmd.format(geno_path+'/'+phe_name, fn, p1, p2,str(500), './clump_result/' + phe_name + '_'+str(500)),))
        ms.append(phe_name)
    s = mp.parallel(mp.run, cmds, num_threads)
    if sum(s) != 0:
        print(','.join(list(np.array(ms)[s]))+' do not  successfully generated clumped file.')
    return s

def merge_qtl_phe(qtl):
    qtl = qtl.sort_values(by=['CHR','qtl_start'])
    merged_phe_qtl = list()
    for index,row in qtl.iterrows():
        if not merged_phe_qtl:
            merged_phe_qtl.append(row)
        else:
            if row['CHR'] != merged_phe_qtl[-1]['CHR']:
                merged_phe_qtl.append(row)
            else:
                if row['qtl_start'] < merged_phe_qtl[-1]['qtl_end'] + 3000000:
                    if row['P'] < merged_phe_qtl[-1]['P']:
                        merged_phe_qtl[-1]['P'] = row['P']
                        merged_phe_qtl[-1]['SNP'] = row['SNP']
                    merged_phe_qtl[-1]['qtl_start'] = min(merged_phe_qtl[-1]['qtl_start'],row['qtl_start'])
                    merged_phe_qtl[-1]['qtl_end'] = max(merged_phe_qtl[-1]['qtl_end'], row['qtl_end'])
                    merged_phe_qtl[-1]['SP2_num'] += row['SP2_num']
                else:
                    merged_phe_qtl.append(row)
    merged_phe_qtl = pd.DataFrame(merged_phe_qtl)
    return merged_phe_qtl

def merge_qtl(qtl):
    qtl = qtl.sort_values(by=['CHR','qtl_start'])
    merged_qtl = pd.DataFrame()
    for index,row in qtl.iterrows():
        if merged_qtl.empty:
            merged_qtl = pd.concat([merged_qtl, row.to_frame().T])
        else:
            qtl_length = row['qtl_end'] - row['qtl_start']
            qtl_ratio = (merged_qtl['qtl_end'] - row['qtl_start']) / qtl_length
            qtl_index = (merged_qtl['CHR'] == row['CHR']) & (qtl_ratio >=0.1)
            if qtl_index.sum()>0:
                peak_dis = (merged_qtl.loc[qtl_index,'SNP'].apply(lambda x: int(x.split('_')[-1])) - int(row['SNP'].split('_')[-1])).abs()
                if (peak_dis <= 2000000).sum()==0:
                    merged_qtl = pd.concat([merged_qtl, row.to_frame().T])
                else:
                    merged_qtl_index = peak_dis[qtl_index].idxmin()
                    if merged_qtl.loc[qtl_index,:].shape[0]>1:
                        print(merged_qtl.loc[qtl_index,:])
                        print(row)
                        print(merged_qtl_index)
                    if merged_qtl.loc[merged_qtl_index,'P'] > row['P']:
                        merged_qtl.loc[merged_qtl_index,'P'] = row['P']
                        merged_qtl.loc[merged_qtl_index,'SNP'] = row['SNP']
                    merged_qtl.loc[merged_qtl_index,'qtl_start'] = min(merged_qtl.loc[merged_qtl_index,'qtl_start'],row['qtl_start'])
                    merged_qtl.loc[merged_qtl_index,'qtl_end'] = max(merged_qtl.loc[merged_qtl_index,'qtl_end'], row['qtl_end'])
                    merged_qtl.loc[merged_qtl_index,'SP2_num'] += row['SP2_num']
                    merged_qtl.loc[merged_qtl_index,'phe_name'] = merged_qtl.loc[merged_qtl_index,'phe_name'] + ',' + row['phe_name']

            else:
                merged_qtl = pd.concat([merged_qtl, row.to_frame().T])
    return merged_qtl

def phe_cluster(phe, phe_labeled, n):
    pca = PCA(n_components=1)
    phe_pc1 = pca.fit_transform(phe)
    phe_corr = phe.corrwith(pd.Series(phe_pc1[:,0],index=phe.index)).abs()
    if (phe_corr >= 0.6).all():
        phe_labeled.loc[phe_corr.index,'label'] = n
        return pd.DataFrame(phe_pc1,index=phe.index,columns=['cluster'+str(n)+'_PC1']), phe_labeled, n+1
    else:
        phe_pc1= pd.DataFrame()
        while not phe.empty:
            phe_corr = phe.corrwith(pd.Series(pca.fit_transform(phe)[:,0],index=phe.index)).abs()
            if (phe_corr < 0.6).sum()==1:
                if phe_corr.shape[0]==2:
                    phe_pc1 = pd.concat([phe_pc1,phe.loc[:,phe_corr.index]],axis=1)
                else:
                    phe_pc1 = pd.concat([phe_pc1,pd.DataFrame(pca.fit_transform(phe.loc[:,phe_corr>=0.6]),index=phe.index,columns=['cluster'+str(n)+'_PC1'])],axis=1)
                    phe_labeled.loc[phe.loc[:,phe_corr>=0.6].columns,'label'] = n
                    n = n + 1
                    phe_pc1 = pd.concat([phe_pc1,phe.loc[:, phe_corr < 0.6]],axis=1)
                phe = pd.DataFrame()
            else:
                if (phe_corr>=0.6).any():
                    if (phe_corr>=0.6).sum()==1:
                        phe_pc1 = pd.concat([phe_pc1,phe.loc[:,phe_corr>=0.6]],axis=1)
                    else:
                        phe_pc1 = pd.concat([phe_pc1,pd.DataFrame(pca.fit_transform(phe.loc[:,phe_corr>=0.6]),index=phe.index,columns=['cluster'+str(n)+'_PC1'])],axis=1)
                        phe_labeled.loc[phe.loc[:,phe_corr>=0.6].columns,'label'] = n
                        n = n + 1
                        print(phe.loc[:,phe_corr>=0.6].corrwith(pd.Series(pca.fit_transform(phe.loc[:,phe_corr>=0.6])[:,0],index=phe.index)))
                    phe = phe.loc[:,phe_corr < 0.6]
                else:
                    phe_pc1 = pd.concat([phe_pc1,phe],axis=1)
                    phe= pd.DataFrame()
                #phe_corr = phe.corrwith(pd.Series(pca.fit_transform(phe)[:,0],index=phe.index)).abs()
    return phe_pc1,phe_labeled,n



def generate_qtl(clump_result_dir):
    qtl_res = list()
    bad_qtl = list()
    for fn in glob.glob(clump_result_dir.strip('/')+'/*clumped'):
        phe_name = '_'.join(fn.split('/')[-1].split('_')[1:-1])
        clump_result = pd.read_csv(fn,sep='\s+')
        clump_result = clump_result.loc[clump_result.SP2!='NONE',:]
        qtl = clump_result[['CHR','BP','SNP','P','SP2']]
        qtl.loc[:,'SP2_num'] = qtl['SP2'].apply(lambda x: len(x.split(',')))
        qtl.loc[:,'log10P'] = -np.log10(qtl['P'])
        if (qtl['SP2_num'] >= 10).sum() > 0:
            qtl['qtl_start'] = qtl['SP2'].apply(lambda x:int(re.findall(r'_(\d+)',x)[0]))
            qtl['qtl_end'] = qtl['SP2'].apply(lambda x:int(re.findall(r'_(\d+)',x)[-1]))
            qtl['phe_name'] = phe_name
            qtl_filter = qtl.loc[qtl.SP2_num>=5,:]
            mer_qtl_filter = merge_qtl_phe(qtl_filter)
            if mer_qtl_filter.shape[0] < 10:
                qtl = qtl.loc[(qtl.SP2_num>=5) & (qtl.log10P >= (qtl.loc[qtl['SP2_num'].idxmax(), 'log10P' ]-6) / 2 + 6),['CHR','qtl_start','qtl_end','SNP','P','SP2_num','phe_name']]
                mer_qtl = merge_qtl_phe(qtl)
                mer_qtl.loc[:,'qtl_length'] = mer_qtl['qtl_end'] - mer_qtl['qtl_start'] + 1
                mer_qtl = mer_qtl.loc[:,['CHR','qtl_start','qtl_end','SNP','P','SP2_num','qtl_length','phe_name']]
                if mer_qtl.shape[0] < 4:
                    qtl_res.append(mer_qtl)
                else:
                    bad_qtl.append(mer_qtl.loc[mer_qtl['SP2_num']>=10,:])
            else:
                bad_qtl.append(mer_qtl_filter.loc[mer_qtl_filter['SP2_num']>=10,['CHR','qtl_start','qtl_end','SNP','P','SP2_num','qtl_length','phe_name']])
    qtl_res = pd.concat(qtl_res)
    qtl_res = qtl_res.loc[qtl_res['SP2_num']>=10,:]
    bad_qtl = pd.concat(bad_qtl)
    return qtl_res, bad_qtl

def phe_PCA(omics_phe, qtl):
    qtl_uniq = pd.DataFrame()
    for phe_name in qtl['phe_name'].unique():
        qtl_sub = qtl.loc[qtl.phe_name==phe_name,:]
        if qtl_sub.shape[0]==1:
            qtl_uniq = pd.concat([qtl_uniq,qtl_sub])
        else:
            qtl_sub = qtl_sub.sort_values(by=['SP2_num'],ascending=False)
            if qtl_sub.iloc[0,:]['SP2_num'] / qtl_sub.iloc[1,:]['SP2_num'] > 2:
                qtl_uniq = pd.concat([qtl_uniq,qtl_sub.iloc[0,:].to_frame().T])
            else:
                qtl_uniq = pd.concat([qtl_uniq,qtl_sub.loc[qtl_sub['P'].idxmin(),:].to_frame().T])
    merged_qtl_uniq = merge_qtl(qtl_uniq)
    merged_qtl_uniq.loc[:,'phe_name_num'] = merged_qtl_uniq['phe_name'].apply(lambda x:len(x.split(',')))
    merged_qtl_uniq = merged_qtl_uniq.sort_values(by='phe_name_num',ascending=False)
    omics_sig_phe = omics_phe.loc[:,[p.replace('m.z','m/z') for p in qtl['phe_name'].unique()]]
    omics_sig_phe_labeled = omics_sig_phe.T
    omics_sig_phe_labeled.loc[:,'label'] = 0
    omics_phe_pc = pd.DataFrame()
    n=1
    for phe_name in merged_qtl_uniq.loc[merged_qtl_uniq['phe_name_num']>=2,'phe_name']:
        omics_phe_sub = omics_sig_phe.loc[:,[p.replace('m.z','m/z') for p in phe_name.split(',')]]
        omics_phe_sub_pc,omics_sig_phe_labeled,n = phe_cluster(omics_phe_sub, omics_sig_phe_labeled, n)
        omics_phe_pc = pd.concat([omics_phe_pc, omics_phe_sub_pc],axis=1)
    clustered_omics_phe = pd.merge(omics_phe_pc.loc[:,~omics_phe_pc.columns.isin(omics_sig_phe.columns)],omics_sig_phe_labeled.loc[omics_sig_phe_labeled.label==0,:].drop('label',axis=1).T,left_index=True,right_index=True)
    return clustered_omics_phe,omics_sig_phe_labeled
