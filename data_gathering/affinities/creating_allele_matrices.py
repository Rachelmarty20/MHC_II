import pandas as pd
import sys


def main(category):

    alleles = [x.strip() for x in open('/cellar/users/ramarty/Data/hla_ii/presentation/other/netMHCIIpan_alleles.txt').readlines() if x.strip()[:2] != 'H-']

    mutations = [x.strip() for x in open('/cellar/users/ramarty/Data/hla_ii/presentation/residues/{0}.txt'.format(category)).readlines()]
    df = pd.DataFrame({'mutation': mutations})

    for i, allele in enumerate(alleles):
        BR = []
        aff = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/affinities/{0}/{1}.csv'.format(category, allele), sep='\t', skiprows=1)
        for mutation in mutations:
            BR.append(aff[aff.ID == mutation].Rank.min())
        df[allele] = BR
    df.to_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.csv'.format(category))



###########################################  Main Method  #####################################

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print 'Wrong number of arguments.'
        sys.exit()
    main(sys.argv[1])
    sys.exit()