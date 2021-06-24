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

echo "export PATH=`pwd`/utils:\$PATH" >> ~/.bashrc
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

echo "export PATH=`pwd`/utils:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

## A toy try
### Downloading example data
```
git clone https://github.com/liusy-jz/MODAS_data.git
```
First check the integrity of the downloaded data, MODAS_data contains five folders, namely  agronomic_traits, genotype, metabolome, transcriptome and example_data. The example folder contains sample data for MODAS, while other folders contain the omics data used in the article.<br/><br/>
Then, enter the MODAS_data directory,
```
cd MODAS_data
```
### Generate pseudo-genotype files
```
MODAS.py genoidx -g example_data/example_geno -genome_cluster -o example_geno
```
pseudo-genotype files generated by genoidx subcommand will be saved as example_geno.genome_cluster.csv.
### Prescreen candidate genomic regions for omics data
The prescreen subcommand uses genome-wide genotype files to calculate the kinship matrix, first extract genotype files by:
```
tar -xvf genotype/chr_HAMP_genotype.tar.gz
```
Then, the pseudo-genotype file example_geno.genome_cluster.csv generated by genoidx and the example_phe.csv file under the example_data folder are used for prescreen analysis,
```
MODAS.py prescreen -g ./genotype/chr_HAMP -genome_cluster example_geno.genome_cluster.csv -phe example_data/example.phe.csv -o example
```
Prescreen subcommand  generates two files including example.sig_omics_phe.csv containing phenotype data and example.phe_sig_qtl.csv containing candidate genomic regions of phenotype.
### Perform regional association analysis to identify QTLs
The prescreen subcommand outputs are used for regional association analysis,
```
MODAS.py regiongwas -g ./genotype/chr_HAMP -phe example.sig_omics_phe.csv -phe_sig_qtl example.phe_sig_qtl.csv -o example
```
Regiongwas subcommand generates two QTL files including example.region_gwas_qtl_res.csv containing reliable QTL results and example.region_gwas_bad_qtl_res.csv containing unreliable QTL results.
### Perform Mendelian randomization analysis
```
MODAS.py mr -g ./genotype/chr_HAMP -exposure ./example_data/test.exp.csv -outcome agronomic_traits/blup_traits_final.new.csv -qtl example_data/example_qtl_res.csv -mlm -o example
```
The results of Mendelian randomization analysis are saved as example.MR.csv.
### MR-based network analysis
MR-based network analysis is carried out by the parameter net of mr subcommand. It uses transcriptome data for subnetwork modules analysis,
```
MODAS.py mr -g ./chr_HAMP/genotype -exposure ./example_data/network_example.exp.csv -outcome ./example_data/network_example.exp.csv -qtl example_data/network_example_qtl.csv -mlm -net -o network_example
```
Network analysis generated four files, including network_example.MR.csv containing gene pairs with MR effect, network_example.edgelist containing gene pairs with weight, network_example.cluster_one.result.csv containing all identified subnetwork modules, network_example.sig.cluster_one.result.csv containing significant subnetwork modules.


## Document
detail in [https://modas-bio.github.io/](https://modas-bio.github.io/)
