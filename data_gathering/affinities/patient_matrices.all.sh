#! /bin/csh
#$ -V
#$ -S /bin/csh
#$ -o /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -e /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -cwd
#$ -t 1-3
#$ -l h_vmem=2G
#$ -tc 50
#$ -l long
set categories=(tsgenes oncogenes indels)
set conditions=(mut mut mut)

set category=$categories[$SGE_TASK_ID]
set condition=$conditions[$SGE_TASK_ID]

date
hostname

python /cellar/users/ramarty/Projects/hla_ii/data_gathering/affinities/creating_patient_matrices.all.py $category alternate $condition

date
