#Rscript multi_trait_signal_manhattan.R gwas_qtl_res.csv background_snps dir_of_assoc_files title_of_figure
library(data.table)
library("CMplot")

options(datatable.fread.datatable=FALSE)
args <- commandArgs(trailingOnly = TRUE)

qtl <- read.csv(args[1],header=T,stringsAsFactors=F)
background_snps_file <- args[2]
dir_of_assoc_files <- args[3]
title_of_figure <- args[4]

# make background signal scatter plot
tmp_bk <- fread(background_snps_file,stringsAsFactors = F,header = T, sep='\t')
tmp_bk$p_wald <- runif(dim(tmp_bk)[1],1e-5,1e-3)

for(p in unique(qtl$phe_name)){
	fn <- paste(dir_of_assoc_files,"/","tmp_",p,'_plink.assoc.txt',sep='')
	g <- fread(fn,stringsAsFactors = F,header = T, sep='\t')
	g <- subset(g,p_wald<=1e-2)
	g[g$p_wald<1e-20,'p_wald'] <- 1e-20
	#sig_g <- subset(g,chr==g[g$p_wald==min(g$p_wald),'chr'] & ps >= g[g$p_wald==min(g$p_wald),'ps']-1000000 & ps <= g[g$p_wald==min(g$p_wald),'ps']+1000000)
	signal <- qtl[qtl$phe_name==p,]
	for(n in 1:nrow(signal)){
		sig_g <- subset(g,chr==signal[n,'CHR'] & ps >= signal[n,'qtl_start'] & ps <= signal[n,'qtl_end'])
        tmp_bk[match(sig_g$rs,tmp_bk$rs),'p_wald'] <- sig_g$p_wald
	}	
	if(is.null(sig_g)){
		next
	}
	#bk[match(sig_g$rs,bk$rs),'p_wald']<-sig_g$p_wald
}

thresholdi <- c(1/dim(tmp_bk)[1], 1e-06, 1e-05)
tmp_bk <- tmp_bk[c('rs','chr','ps','p_wald')]
lim=-log10(min(tmp_bk$p_wald))+2
colnames(tmp_bk) <- c("SNP","Chromosome","Position",title_of_figure)
CMplot(tmp_bk, plot.type=c("m"), col=c("grey30","grey60"), ylim=c(0,lim), threshold=thresholdi,
        cex=c(0.5,0.5,0.5), signal.cex=c(0.5,0.5,0.5),
        threshold.col=c("red",'green','blue'), chr.den.col=NULL, amplify=TRUE, 
        signal.pch=c(19,19,19),
        signal.col=c("red",'green','blue'), multracks=FALSE, LOG10=TRUE,file='jpg')
