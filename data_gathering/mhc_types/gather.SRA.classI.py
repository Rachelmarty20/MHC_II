import pandas as pd
import cPickle as pickle
import samples
import numpy as np

# get TCGA types
samples = list(set([x.strip() for x in open('/cellar/users/ramarty/Data/hla_ii/hla_types/samples.SRA.txt').readlines() if x != '\n']))

# Dictionary with list of alleles for each patient

all_patient_dictionary = {}
for sample in samples:
    print sample

    try:
        patient_dictionary = {}
        directory = '/data/nrnb03/users/ramarty/{0}'.format(sample)

        f = open('{0}/result.txt'.format(directory)).readlines()[1].split('\t')

        patient_dictionary = {}
        for g, i in zip(['A', 'B', 'C'], [1, 3, 5]):
            patient_dictionary[g+'_allele1'] = '{0}_{1}{2}'.format(f[i][0], f[i][2:4], f[i][5:7])
            patient_dictionary[g+'_allele2'] = '{0}_{1}{2}'.format(f[i+1][0], f[i+1][2:4], f[i+1][5:7])

        all_patient_dictionary[sample] = patient_dictionary
    except:
        print "Patient not typed."

df = pd.DataFrame(all_patient_dictionary).transpose()

df.to_csv('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.SRA.classI.csv')
pickle.dump(all_patient_dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.SRA.classI.p', 'wb'))

# A dictionary with a dictionary for each patient containing the specific alleles

patient_types_df = df
patients = list(patient_types_df.index)

dictionary = {}
for patient in patients:
    patient_dictionary = {}
    try:
        # A
        alleles = list(patient_types_df.ix[patient][['A_allele1', 'A_allele2']])
        patient_dictionary['A'] = alleles
        alleles = list(patient_types_df.ix[patient][['B_allele1', 'B_allele2']])
        patient_dictionary['B'] = alleles
        alleles = list(patient_types_df.ix[patient][['C_allele1', 'C_allele2']])
        patient_dictionary['C'] = alleles

        dictionary[patient] = patient_dictionary
    except:
        # weeds out the patients with the failed typing
        print patient

pickle.dump(dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/SRA.HLA_classI.p', 'w'))

# add something for CLASSI: ~/Data/hla_ii/hla_types/optitype_output.txt