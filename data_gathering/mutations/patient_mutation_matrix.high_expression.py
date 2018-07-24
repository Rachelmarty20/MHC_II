import pandas as pd
import sys

peptide_type = sys.argv[1]

patient_mutations = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.cancer.TCGA.inclusive.mut.csv',
                               index_col=0)
mutations = list(patient_mutations.columns)
samples_with_mutations = list(patient_mutations.index)


samples_combined = pd.read_csv('/cellar/users/ramarty/Data/hla/mutations/processed_mutation_files.expression_greater_25q.csv',
                               index_col=0)

sample_mutation_df = {}
for i, sample in enumerate(samples_with_mutations):
    if i % 100 == 0:
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

matrix.to_csv('/cellar/users/ramarty/Data/hla_ii/presentation/patient_matrices/cancer.patient_mutation.expression_greater_25q.csv')

# copied to /cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.low_allelic_fraction.TCGA.inclusive.mut.csv