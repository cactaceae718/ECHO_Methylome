args<-commandArgs(T)
in_dir <- args[1]
out_dir <- args[2]

library("methylKit")

in_dir="/gpfs/data/proteomics/projects/TRASANDE_LAB/methylcall_percent_fullreport/"
out_dir="/gpfs/data/proteomics/projects/TRASANDE_LAB/methy_call_full_report_flatfile_fn/"

myobj=readRDS(file=paste0(out_dir, "myobj.rds"))

files=list.files(path=in_dir, pattern="methycall_percent_report.txt", recursive = TRUE)
files_list <- list()
for (i in 1:length(files)){
  files_list[[i]] <- files[i]
}
print(files)

id=as.character(sapply(files, function(x) list(strsplit(x,'_')[[1]][1])))
id_list <- list()
for (i in 1:length(id)){
  id_list[[i]] <- id[i]
}

#Create QC Report#
for (i in c(1:length(id_list))) {
  pdf(file=paste0(out_dir,id_list[i],"MethylationStats.pdf"),width=14,height=8) ; getMethylationStats(myobj[[i]],plot=T,both.strands=T); dev.off()
  pdf(file=paste0(out_dir,id_list[i],"CpGcoverageStats.pdf"),width=14,height=8) ; getCoverageStats(myobj[[i]],plot=T,both.strands=T); dev.off()
  }
#Create QC Report#                         