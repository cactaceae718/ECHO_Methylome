---
title: "methyl_n169_treatment_toba"
author: "YC"
date: "6/8/2019"
output: html_document
---

in_dir <- "/gpfs/data/proteomics/projects/TRASANDE_LAB/methylcall_percent_fullrepor"
out_dir <- "/gpfs/data/proteomics/projects/TRASANDE_LAB/methycall_full_report_flatfile_fn/"

dir.create(out_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

files = list.files(path=in_dir, pattern="methycall_percent_report.txt",recursive = TRUE)
print(files)
files_list <- list()
for (i in 1:length(files)){
  files_list[[i]] <- files[i]
}
print(files_list)

id <- as.character(sapply(files, function(x) list(strsplit(x,'_')[[1]][1])))
id_list <- list()
for (i in 1:length(id)){
  id_list[[i]] <- id[i]
}
print(id_list)

#get nonsmoker and smoker ID
ds_Mdna <- read.table('ds_Mdna.txt', header = T, sep = "\t", row.names = 1)
id_smoker <- na.omit(rownames(ds_Mdna)[ds_Mdna$tobacco_binary=="Yes"])
id_nonsmoker <- na.omit(rownames(ds_Mdna)[ds_n461$tobacco_binary=="No"])

library("methylKit")
#Create MethyKit Object by keeping all loci mincov=0. default min.cov=10#
myobj <- methRead(files_list, sample.id=id_list, assembly="hg38", treatment=rep(1,length(id_list)), context="CpG",dbtype = "tabix", dbdir = getwd(), mincov=0) 
saveRDS(myobj,file=paste0(out_dir,"myobj.rds"))
#myobj <- readRDS(paste0(out_dir, "myobj.rds"))

#re-organize samples and treatment vector within methylRawListDB objects
myobj_toba=reorganize(myobj, sample.ids=c(id_toba, id_nontoba), treatment=c(rep(0, length(id_toba)), rep(1, length(id_nontoba))))

#N=169, full model(treatment toba vs non_toba), unite(all methyl profile by keeping regions/bases that are covered in all samples)
meth_toba=unite(myobj_toba, mc.cores=2)
saveRDS(meth_toba, file = "meth_toba.rds")

#Extract %methylation values from methylBase object
perc_meth <- percMethylation(meth_toba)
#write.table(perc_meth, file = "perc_meth.txt", row.names = F, col.names = T, sep = "\t")

#methylation differentional
myDiff_toba=calculateDiffMeth(meth_toba)

###### select the bases that have q-value<0.05 and percent methylation difference > 15% ######
diff15p <- getMethylDiff(myDiff_toba,difference=15,qvalue=0.05)
library("genomation")
df_diff15p <- as.data.frame(as(diff15p, "GRanges"))
#re-format as input for plot(qqman pkg)
df_diff15p$snpID <- paste(df_diff15p$seqnames,"-",df_diff15p$start, sep = "")
df_diff15p_rev <- df_diff15p[,c(9,1,2,7,8)] #remove pseudo chromosome and re-order column for plotting
colnames(df_diff15p_rev) <- c("SNP", "CHR", "BP", "P", "meth.diff")
#re-format to numeric
df_diff15p_rev$CHR <- as.integer(substr(as.character(df_diff15p_rev$CHR), 4,5))
#DMS==Differential Methylated Sites
DMS <-df_diff15p_rev$SNP

library("CMplot") # perform Manhattan Plot
#chromosomal methylation coverage p=0.05, meth.diff=15
CMplot(df_diff15p_rev[,1:4], plot.type=c("c","d"), r=1.6, cir.legend=TRUE,
        outward=TRUE, cir.legend.col="black", cir.chr.h=1.1 ,chr.den.col="orange", file="jpg",
        memo="", dpi=300, chr.labels=paste("Chr",c(1:22), sep=""))

###### all the bases that have q-value<0.05 ######
#re-format as input for plot(qqman pkg)
diffq=getMethylDiff(myDiff_toba, difference=0, qvalue=0.05)
df_diffq<- as.data.frame(as(diffq, "GRanges"))
df_diffq$snpID <- paste(df_diffq$seqnames,"-",df_diffq$start, sep = "")
df_diffq_rev <- df_diffq[c(4:351,354:601),c(9,1,2,7,8)] #remove pseudo chromosome and re-order column for plotting
colnames(df_diffq_rev) <- c("SNP", "CHR", "BP", "P", "meth.diff")
df_diffq_rev$CHR <- as.integer(substr(as.character(df_diffq_rev$CHR), 4,5))

# chromosomal methylation coverage p=0.05, meth.diff=0 #
CMplot(df_diffq_rev[,1:4], plot.type=c("c","d"), r=1.6, cir.legend=TRUE,
        outward=TRUE, cir.legend.col="black", cir.chr.h=1.1 ,chr.den.col="orange", file="jpg",
        memo="", dpi=300, chr.labels=paste("Chr",c(1:22), sep=""), threshold = 0.01)

library("qqman")
#cutoff p-value<=0.05, Bonferroni correction=p-value/# of observations
manhattan(df_diffq_rev, main = "Manhattan Plot", logp = T, suggestiveline = -log10(5e-2), genomewideline = -log10(0.05/92700), highlight = DMS, annotatePval = 0.01, annotateTop = T)
qq(df_diff15p$P, main = "Q-Q plot of diff_methy p-values")

###### all the bases that have aligned ######
library(genomation)
df_myDiff <- as.data.frame(as(myDiff_toba,"GRanges"))
df_myDiff$snpID <- paste(df_myDiff$seqnames,"-",df_myDiff$start, sep = "")
df_myDiff_rev <- df_myDiff[,c(9,1,2,7,8)] #remove pseudo chromosome and re-order column for plotting
colnames(df_myDiff_rev) <- c("SNP", "CHR", "BP", "P", "meth.diff")
#re-format to satisfy library(qqman)
df_myDiff_rev$CHR <- as.integer(substr(as.character(df_myDiff_rev$CHR), 4,5))

library("CMplot")
CMplot(df_myDiff_rev[1:89438,1:4], plot.type=c("c","d"), r=1.6, cir.legend=TRUE,
        outward=TRUE, cir.legend.col="black", cir.chr.h=1.1 ,chr.den.col="orange", file="jpg",
        memo="", dpi=300, chr.labels=paste("Chr",c(1:22), sep=""), threshold = 0.01)

library("qqman")
manhattan(df_myDiff_rev[1:89438,], main = "Manhattan Plot", suggestiveline = -log10(5e-2), genomewideline = -log10(0.05/92700), highlight = DMS, annotatePval = 0.01, annotateTop = T, col = c("pink2", "orange"))
#only plot chrosome 4 and 8 
manhattan(subset(df_myDiff_rev, CHR == 8 | CHR == 4, highlight = DMS))
#perform qq plot
qq(df_myDiff_rev$P, main = "Q-Q plot of diff_methy p-values")

===================================================================================================
###### Annotation based on limma 
library(genomation)
#read the gene BED file
gene_obj <- readTranscriptFeatures("/Users/cactaceae/Desktop/Methylomics/refseq.hg38.bed.txt")

#annotate differentially methylated CpGs with promoter/exon/intron using annotation data
diffAnn <- annotateWithGeneParts(as(Diff15p,"GRanges"),gene_obj)

#diffAnn_overlap <- annotateWithGeneParts(as(Diff25p,"GRanges"),gene_obj,intersect.chr=TRUE)
getAssociationWithTSS(diffAnn)
getTargetAnnotationStats(diffAnn, percentage=TRUE, precedence=TRUE)
#diffAnn@precedence == TargetAnnotationStats

plotTargetAnnotation(diffAnn,precedence=TRUE, main="differential methylation annotation") #TRUE percentage of annotation feature. FALSE, number of annotation features
getFeatsWithTargetsStats(diffAnn,percentage=F)

#diffAnn@perc.of.OlapFeat==getFeatsWithTargetsStats
                          
#read the shores and flanking regions and name the flanks as shores and CpG islands as CpGi
cpgi_obj=readFeatureFlank("/Users/cactaceae/Desktop/Methylomics/cpgIslandExt_rev.txt",feature.flank.name=c("CpGi","shores"))
diffCpGann_10p=annotateWithFeatureFlank(as(Diff10p,"GRanges"), cpgi_obj$CpGi,cpgi_obj$shores, feature.name="CpGi",flank.name="shores")
plotTargetAnnotation(diffCpGann_10p,col=c("green","gray","white"), main="differential methylation annotation")                          
```
