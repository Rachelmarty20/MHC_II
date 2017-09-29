#! /bin/csh
#$ -V
#$ -S /bin/csh
#$ -o /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -e /cellar/users/ramarty/Data/hla_ii/sge-system_files
#$ -cwd
#$ -t 1-195
#$ -l h_vmem=1G
#$ -tc 200
#$ -l long
set alleles=(DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101 DRB1_0401 DRB1_1301 DRB1_0401 DRB1_1301 DRB1_0301 DRB1_1101 DRB1_0901 DRB1_1001 DRB1_0701 DRB1_1501 DRB5_0101 DRB1_0101 DRB1_0701 DRB1_0101 DRB1_1101)
set donors=(DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE DonorB DonorB DonorC DonorC DonorA DonorA DonorF DonorF DonorG DonorG DonorG DonorD DonorD DonorE DonorE)
set lengths=(13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25)

set allele=$alleles[$SGE_TASK_ID]
set donor=$donors[$SGE_TASK_ID]
set length=$lengths[$SGE_TASK_ID]

date
hostname

echo starting netMHCIIpan
/cellar/users/ramarty/programs/netMHCIIpan-3.1/netMHCIIpan -f /cellar/users/ramarty/Data/hla_ii/validation/ciudad/fasta_files/$donor.$length.fa -a $allele -length $length -xls -xlsfile /cellar/users/ramarty/Data/hla_ii/validation/ciudad/affinities/$donor.$allele.$length.csv > /cellar/users/ramarty/Data/hla_ii/sge-system_files/$donor.$allele.$length.txt
/cellar/users/ramarty/programs/netMHCIIpan-3.1/netMHCIIpan -f /cellar/users/ramarty/Data/hla_ii/validation/ciudad/fasta_files/$donor.random.$length.fa -a $allele -length $length -xls -xlsfile /cellar/users/ramarty/Data/hla_ii/validation/ciudad/affinities/$donor.$allele.random.$length.csv > /cellar/users/ramarty/Data/hla_ii/sge-system_files/$donor.$allele.random.$length.txt

date
