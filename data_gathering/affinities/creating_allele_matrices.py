import pandas as pd
import sys

category = sys.argv[1]

alleles = [x.split('\t')[0].strip().split('*')[0] + '_' + ''.join(x.split('\t')[0].strip().split('*')[1].split(':')) for x in open('/cellar/users/ramarty/Data/hla_ii/presentation/other/alleles.txt').readlines()[1:]]

mutations = [x.strip() for x in open('/cellar/users/ramarty/Data/hla_ii/presentation/residues/{0}.txt'.format(category)).readlines()]
df = pd.DataFrame({'mutation': mutations})

for i, allele in enumerate(alleles):
    if i % 100 == 0:
        print i
    BR = []
    aff = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/affinities/{0}/{1}.csv'.format(category, allele), sep='\t', skiprows=1)
    for mutation in mutations:
        BR.append(aff[aff.ID == mutation].Rank.min())
    df[allele] = BR
df.to_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.csv'.format(category))