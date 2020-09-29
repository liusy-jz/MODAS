import pandas as pd
import numpy as np
import modas.multiprocess as mp
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


def merge_qtl(qtl):
    qtl = qtl.sort_values(by=['CHR','BP'])
    merged_qtl = list()
    for index,row in qtl.iterrows():
        if not merged_qtl:
            merged_qtl.append(row)
        else:
            if row['CHR'] != merged_qtl[-1]['CHR']:
                merged_qtl.append(row)
            else:
                if row['BP'] - merged_qtl[-1]['BP'] <= 1000000:
                    if row['P'] < merged_qtl[-1]['P']:
                        merged_qtl[-1]['P'] = row['P']
                        merged_qtl[-1]['BP'] = row['BP']
                        merged_qtl[-1]['SNP'] = row['SNP']
                    merged_qtl[-1]['SP2_num'] += row['SP2_num']
                    merged_qtl[-1]['SP2']+= ',' + row['SP2']
                else:
                    merged_qtl.append(row)
    merged_qtl = pd.DataFrame(merged_qtl)
    return merged_qtl


def generate_qtl(clump_result_dir):
    qtl_res = list()
    for fn in glob.glob(clump_result_dir.strip('/')+'/*clumped'):
        phe_name = '_'.join(fn.split('/')[-1].split('_')[1:-1])
        clump_result = pd.read_csv(fn,sep='\s+')
        clump_result = clump_result.loc[clump_result.SP2!='NONE',:]
        qtl = clump_result[['CHR','BP','SNP','P','SP2']]
        qtl['SP2_num'] = qtl['SP2'].apply(lambda x: len(x.split(',')))
        qtl['qtl_start'] = qtl['SP2'].apply(lambda x:int(re.findall(r'_(\d+)',x)[0]))
        qtl['qtl_end'] = qtl['SP2'].apply(lambda x:int(re.findall(r'_(\d+)',x)[-1]))
        qtl['phe_name'] = phe_name
        qtl['qtl_length'] = qtl['qtl_end'] - qtl['qtl_start'] + 1
        qtl = qtl.loc[qtl.SP2_num>=10,['CHR','qtl_start','qtl_end','SNP','P','qtl_length','phe_name']]
        qtl_res.append(qtl)
    qtl_res = pd.concat(qtl_res)
    return qtl_res
