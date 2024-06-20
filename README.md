# MODAS
## Introduction
MODAS (Multi-Omics Data Association Study toolkit) is an efficient software for high-dimensional omics data association analysis, featuring five main characteristics.

+ MODAS employs a novel whole-genome association study (GWAS) strategy to handle high-dimensional omics data. Specifically, MODAS generates pseudo-genotype files to reduce the dimensionality of the original SNP-based genotype data. It then performs block-based GWAS to identify significantly associated genomic regions (SAGR), and finally conducts SNP-based association analysis on the SAGR to obtain QTLs for high-dimensional omics data.

+ MODAS introduces Mendelian randomization algorithms to infer genetic regulatory relationships between molecular QTLs and phenotypic traits. This method helps to establish biological hypotheses and uncover potential genetic mechanisms.

+ MODAS uses contrastive PCA algorithm to uncover differential components in population omics data between different conditions. This functionality is particularly useful for identifying agronomically important genes involved in condition comparisons, such as genes related to crop stress tolerance.

+ MODAS employs an image matching algorithm to integrate multi-omics molecular QTLs, aiding in the elucidation of the genetic mechanisms underlying complex traits.

+ MODAS provides an HTML-based web interface for visualizing and querying GWAS results from omics data, facilitating further gene function exploration.
## Installation
### Install via pip
```
#Depends: R(>=3.6), python(>=3.8)
pip install MODAS
```
### Installation using conda
```
conda create -n modas python=3.8 -y
conda install -y -c conda-forge r-rcppeigen r=3.6 rpy2
$CONDA_PREFIX/bin/pip install MODAS
```
## New features in MODAS2
### Identifying stress-responsive Molecular QTL through contrastive PCA and Two-Way ANOVA
MODAS2 innovatively employs the contrastive PCA algorithm to separate the stress-response effects of natural variation on molecular traits from the genetic effects of natural variation on molecular traits. The separated stress response effects are then used as stress response indices for identifying stress-responsive molecular QTL. After identifying QTL using the stress response indices, the significance of the stress response for these QTL is further assessed through two-way ANOVA, resulting in the final stress-responsive molecular QTL. MODAS2 uses the subcommand `contrast` to identify stress-responsive molecular QTL. The command line is as follows:
```
MODAS contrast -stress_phe example_data/test_salt.phe.csv -control_phe example_data/test_control.phe.csv -gwas -g example_data/example_geno.contrast -genome_cluster example_data/example_contrast.genome_cluster.csv -p 10 -o example
```
`contrast` subcommand generates five files including `example.scpca_pc.phe.csv`, `example.scpca_pc.beta_test.csv`, `example.scpca_pc.normalized_phe.csv`, `example.region_gwas_qtl_res.anova.csv`, `example.region_gwas_qtl_res.anova.sig.csv` and `example.region_gwas_bad_qtl_res.csv`, which contain the stress response indices of molecular traits, statistical test results of the stress response indices, normalized stress response indices of molecular traits, stress response indices QTL results including two-way ANOVA P values, stress-responsive molecular QTL results and unreliable molecular QTL results, respectively.
### Multi-Trait QTL Colocalization Analysis Based on Image Matching Algorithms
QTL colocalization analysis is an effective method for integrating functional information across different traits. However, existing methods rarely perform multi-trait QTL colocalization analysis. MODAS2 employs an image matching algorithm to score the colocalization degree between pairs of QTLs, and then uses clustering algorithms to quickly achieve multi-trait QTL colocalization analysis. Multi-trait QTL colocalization analysis can be performed using the `coloc` subcommand of MODAS2. The command line is as follows:
```
MODAS coloc -qtl example_data/test_coloc.qtl.csv -g example_data/example_geno.contrast -gwas_dir example_data/gwas_coloc_test/ -p 6 -o example
```
`coloc` subcommand generates three files including `example.coloc_res.csv`, `example.coloc_pairwise.csv` and `example.dis_res.csv`, which contain the results of the co-localized QTL clusters, the pairwise QTL co-localization results, and the similarity results between pairwise QTLs, respectively.

__Note__: Sample data for MODAS2 can be downloaded via DOI [https://zenodo.org/doi/10.5281/zenodo.11951520](https://zenodo.org/doi/10.5281/zenodo.11951520).
## MODAS analytical pipeline
### Downloading example data
`MODAS_data` containing sample data for MODAS and omics data used in the article uploaded by Git extension Git Large File Storage (LFS), first download Git LFS from [https://git-lfs.github.com/](https://git-lfs.github.com/),  and place the `git-lfs` binary on your system’s executable `$PATH` or equivalent, then set up Git LFS for your user account by running:
```
git lfs install
```
next download `MODAS_data` by running:
```
git clone https://github.com/liusy-jz/MODAS_data.git
```
When the download is complete, first check the integrity of the downloaded data, `MODAS_data` contains five folders, namely `agronomic_traits`, `genotype`, `metabolome`, `transcriptome` and `example_data`, also contains a gene annotaion file for maize. The example folder contains sample data for MODAS, while other folders contain the omics data used in the article.<br/><br/>
Then, enter the `MODAS_data` directory,
```
cd MODAS_data
```
### Generate pseudo-genotype files
```
MODAS genoidx -g example_data/example_geno -genome_cluster -o example_geno
```
Pseudo-genotype files generated by `genoidx` subcommand will be saved as `example_geno.genome_cluster.csv`.
### Prescreen candidate genomic regions for omics data
The `prescreen` subcommand uses genome-wide genotype files to calculate the kinship matrix, first extract genotype files by:
```
tar -xvf genotype/chr_HAMP_genotype.tar.gz
```
Then, the pseudo-genotype file `example_geno.genome_cluster.csv` generated by `genoidx` and the `example_phe.csv` file under the `example_data` folder are used for prescreen analysis,
```
MODAS prescreen -g ./chr_HAMP -genome_cluster example_geno.genome_cluster.csv -phe example_data/example.phe.csv -o example
```
`prescreen` subcommand generates two files including `example.sig_omics_phe.csv` containing phenotype data and `example.phe_sig_qtl.csv` containing candidate genomic regions of phenotype.
### Perform regional association analysis to identify QTLs
The `prescreen` subcommand outputs are used for regional association analysis,
```
MODAS regiongwas -g ./chr_HAMP -phe example.sig_omics_phe.csv -phe_sig_qtl example.phe_sig_qtl.csv -o example
```
`regiongwas` subcommand generates two QTL files including `example.region_gwas_qtl_res.csv` containing reliable QTL results and `example.region_gwas_bad_qtl_res.csv` containing unreliable QTL results.
### Perform Mendelian randomization analysis
```
MODAS mr -g ./chr_HAMP -exposure ./example_data/example.exp.csv -outcome agronomic_traits/blup_traits_final.new.csv -qtl example_data/example_qtl_res.csv -mlm -o example
```
The results of Mendelian randomization analysis are saved as `example.MR.csv`.
### MR-based network analysis
MR-based network analysis is carried out by the parameter `-net` of `mr` subcommand. It uses transcriptome data for subnetwork modules analysis,
```
MODAS mr -g ./chr_HAMP -exposure ./example_data/network_example.exp.csv -outcome ./example_data/network_example.exp.csv -qtl example_data/network_example_qtl.csv -mlm -net -o network_example
```
Network analysis generated four files, including `network_example.MR.csv` containing gene pairs with MR effect, `network_example.edgelist` containing gene pairs with weight, `network_example.cluster_one.result.csv` containing all identified subnetwork modules, `network_example.sig.cluster_one.result.csv` containing significant subnetwork modules.
### co-associated gene analysis
Co-associated genes analysis is not a modas function. It is implemented by script `co-associated.py`. The analysis command line is as follows:
```
python3 example_data/co-associated.py example_data/co_associated.test.pvalue.csv co-associated_test
```
Then, a file containing co-associated gene labels and a heatmap showing relationship between co-associated genes are saved as `co-associated_test.cluster.csv` and `co-associated_test.cluster.heatmap.pdf`.
## Document
detail in [https://modas-bio.github.io/](https://modas-bio.github.io/)
