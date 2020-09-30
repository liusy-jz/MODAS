# MODAS
MODAS: Multi-omics data association study

## installation
```
conda create -n modas python=3.7 -y
conda activate modas
python setup.py build
python setup.py install

pip install cffi    # if the version of cffi less than required version
conda install -y -c conda-forge r-rcppeigen r=3.6
Rscript -e 'install.packages("Matrix")'
Rscript -e 'install.packages("bigsnpr")'
```

## A toy try
```
MODAS.py genoidx -g ./RNA_seq_162 -clump -genome_cluster -p 10 -o RNA_seq_162
MODAS.py prescreen -g ./RNA_seq_162_clump -genome_cluster RNA_seq_162.genome_cluster.csv -phe salt_RIL_162.filter.boxcox.normalized_phe.csv -lmm_suggest_pvalue 1e-5 -p 20 -o salt_RNA_seq_162_boxcox
MODAS.py localgwas -g ./RNA_seq_162_clump -phe salt_RNA_seq_162_boxcox.sig_omics_phe.csv -phe_sig_qtl salt_RNA_seq_162_boxcox.phe_sig_qtl.csv -p1 1e-7 -p2 1e-6 -p 20 -o salt_RNA_seq_162_boxcox
```
