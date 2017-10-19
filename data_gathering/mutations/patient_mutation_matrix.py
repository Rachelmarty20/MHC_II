import pandas as pd
import sys

peptide_type = sys.argv[1]

mutations = [x.strip() for x in open('/cellar/users/ramarty/Data/hla_ii/presentation/residues/{0}.txt'.format(peptide_type))]

samples_combined = pd.read_csv('/cellar/users/ramarty/Data/hla/mutations/processed_mutation_files.full_tcga.csv', index_col=0)
samples_with_mutations = list(samples_combined['Barcode'].unique())
samples_combined = samples_combined[(samples_combined['Barcode'].isin(samples_with_mutations)) & (samples_combined.combined.isin(mutations))]

sample_mutation_df = {}
for i, sample in enumerate(samples_with_mutations):
    print i
    mutation_array = []
    for mutation in mutations:
        if len(samples_combined[(samples_combined['Barcode'] == sample) & (samples_combined['combined'] == mutation)]) > 0:
            mutation_array.append(1)
        else:
            mutation_array.append(0)
    sample_mutation_df[sample] = mutation_array

matrix = pd.DataFrame(sample_mutation_df)
matrix.index = mutations
matrix = matrix.transpose()

matrix.to_csv('/cellar/users/ramarty/Data/hla_ii/presentation/patient_matrices/{0}.patient_mutation.csv'.format(peptide_type))
