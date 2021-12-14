if(!require(gplots)) install.packages("gplots")
if(!require(LDheatmap)) install.packages('https://cran.r-project.org/src/contrib/Archive/LDheatmap/LDheatmap_0.99-7.tar.gz',repos=NULL,type='source')
if(!require(genetics)) install.packages("genetics")
if(!require(ape)) install.packages("ape")
if(!require(compiler)) install.packages("compiler")

if(!require(EMMREML)) install.packages("EMMREML")
if(!require(scatterplot3d)) install.packages("scatterplot3d")

if(!'multtest'%in% installed.packages()[,"Package"]){
    if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
    BiocManager::install(c("multtest", 'snpStats'))
}
