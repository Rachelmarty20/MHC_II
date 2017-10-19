#! /bin/csh
#$ -V
#$ -S /bin/csh
#$ -o /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -e /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -cwd
#$ -t 1-1
#$ -l h_vmem=2G
#$ -tc 50
#$ -l long
set categories=(germline)

set category=$categories[$SGE_TASK_ID]

date
hostname

python /cellar/users/ramarty/Projects/hla_ii/data_gathering/affinities/creating_allele_matrices.py $category

date
