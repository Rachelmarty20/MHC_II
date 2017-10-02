#! /bin/csh
#$ -V
#$ -S /bin/csh
#$ -o /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -e /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -cwd
#$ -t 1-13
#$ -l h_vmem=1G
#$ -tc 50
#$ -l long
set lengths=(13 14 15 16 17 18 19 20 21 22 23 24 25)

set length=$lengths[$SGE_TASK_ID]

date
hostname

python /cellar/users/ramarty/Projects/hla_ii/data_gathering/affinities/validation/process_ms_peptides.py $length

date
