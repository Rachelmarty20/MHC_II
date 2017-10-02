#! /bin/csh
#$ -V
#$ -S /bin/csh
#$ -o /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -e /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -cwd
#$ -t 1-4
#$ -l h_vmem=8G
#$ -tc 50
#$ -l long
set categories=(DR DP DQ all)

set category=$categories[$SGE_TASK_ID]

date
hostname

Rscript --vanilla /cellar/users/ramarty/Projects/hla_ii/data_analysis/patient_selection/model.R $category > ~/Data/hla_ii/generated_data/output.$category.txt

date
