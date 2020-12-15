# This folder has Shell and R scripts for HPC NYU Langone Phoenix (SGE)
# Project: ECHO_EWAS Methylome Data
This script is for targeted bisulfite sequencing (Methyl-Seq) which is a next-generation sequencing technology used to provide the status of single-base resolution of methylated C by treating the DNA with sodium bisulfite before sequencing. 

Analysis workflow
1.	Quality control analysis: FASTQC
2.	Alignment: DRAGEN Methylation Pipeline 3.2.8 which performs alignment, methyl calling, and calculates alignment and methylation metrics base on BisMark 
3.	Annotation and differential methylation analysis: MethylKit (hg38) (refer Methylkit.r) 
4.	Identification of differentially methylated loci/regions (refer Methylkit.r)
