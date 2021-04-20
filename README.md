# MODAS
MODAS: Multi-Omics Data Association Study toolkit

## installation
### Installation using conda
```
git clone https://github.com/liusy-jz/MODAS.git
cd MODAS
conda create -n modas python=3.7 -y
conda activate modas
python setup.py build
python setup.py install

conda install -y -c conda-forge r-rcppeigen r=3.6 rpy2
Rscript -e 'install.packages(c("data.table", "ggplot2", "ggsignif"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("bigsnpr", dependence=T, repos="https://cloud.r-project.org")'

MODAS_PATH=`pwd`
vi ~/.bashrc
export PATH="$MODAS_PATH/utils:$PATH"
source ~/.bashrc
```

### Mannual Installation
```
#Depends: R(>=3.6), python(>=3.7)

git clone https://github.com/liusy-jz/MODAS.git
cd MODAS
python setup.py build
python setup.py install

pip3 install rpy2
Rscript -e 'install.packages(c("data.table", "ggplot2", "ggsignif"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("bigsnpr",dependence=T, repos="https://cloud.r-project.org")'

MODAS_PATH=`pwd`
vi ~/.bashrc
export PATH="$MODAS_PATH/utils:$PATH"
source ~/.bashrc
```

## A toy try
```
MODAS.py genoidx -g ./chr_HAMP -genome_cluster -p 10 -o chr_HAMP
MODAS.py prescreen -g ./chr_HAMP -genome_cluster ./chr_HAMP.genome_cluster.csv -phe ./E3_log2.normalized_phe.csv -p 20 -o E3_log2
MODAS.py regiongwas -g ./chr_HAMP -phe ./E3_log2.sig_omics_phe.csv -phe_sig_qtl ./E3_log2.phe_sig_qtl.csv -p 20 -o E3_log2
```
