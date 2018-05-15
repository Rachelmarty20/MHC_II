import pandas as pd
import sys


def main(category, gene):

    print category, gene
    alleles = [x.strip() for x in open('/cellar/users/ramarty/Data/hla_ii/presentation/other/netMHCIIpan_alleles.txt').readlines() if x.strip()[:2] != 'H-']

    gene_matrix = []
    for allele in alleles:
        aff = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/whole_gene/affinities/{0}/{1}.csv'.format(category, allele), sep='\t', skiprows=1)
        tmp_aff = aff[aff.ID == gene]
        gene_ranks = []
        for i in range(1, len(tmp_aff) + 15):
            if i <= 15:
                start = 0
                finish = i
            else:
                start = i - 15
                finish = i
            gene_ranks.append(tmp_aff[start:finish].Rank.min())
        gene_matrix.append(gene_ranks)


    df = pd.DataFrame(gene_matrix).transpose()
    df.columns = alleles
    df.to_csv('/cellar/users/ramarty/Data/hla_ii/presentation/whole_gene/allele_matrices/{0}.{1}.csv'.format(category, gene))

###########################################  Main Method  #####################################

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print 'Wrong number of arguments.'
        sys.exit()
    main(sys.argv[1], sys.argv[2])
    sys.exit()