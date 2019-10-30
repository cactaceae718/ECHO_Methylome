#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00
#SBATCH --job-name=methy_cal
#SBATCH --output=methy_call.out

for i in $(ls /gpfs/data/sequence/results/external/NYU/trasandelab/2019-04-15/Methylation_Results)
do
input_dir="/gpfs/data/sequence/results/external/NYU/trasandelab/2019-04-15/Methylation_Results"
sample=$(echo $i | sed 's/_ds.*//')
id=$(echo $sample | cut -d'-' -f 1)
report=$sample.CX_report.txt.gz
output_dir="/gpfs/data/proteomics/projects/TRASANDE_LAB/methylcall_percent_fullreport/"
zcat $input_dir/$i/$report > $output_dir/${id}_methycall_percent_report.txt
done