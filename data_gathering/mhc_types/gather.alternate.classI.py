import pandas as pd
import cPickle as pickle
import samples

# get TCGA types
population_dictionary_stripped = samples.get_population_dictionary()


all_patient_dictionary = {}
for population in population_dictionary_stripped.keys():

    for i, sample in enumerate(population_dictionary_stripped[population]):
        print sample

        patient_dictionary = {}
        directory = '/nrnb/users/ramarty/alternate_pops/{0}/{1}/hlaHD/{2}/result/'.format(population, sample, sample)


        lines = [x.split('\t') for x in open('{0}{1}_final.result.txt'.format(directory, sample)).readlines()[:8]]

        patient_dictionary = {}
        for line in lines:
            print line
            if line[0] in ['A', 'B', 'C']:
                if line[1] == 'Not typed':
                    alleles = ['-', '-']
                elif line[2].strip() == '-':
                    x = line[1]
                    alleles = [line[0]+'_'+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1]]*2
                else:
                    alleles = [line[0]+'_'+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1] for x in line[1:3]]
                patient_dictionary[line[0]+'_allele1'] = alleles[0]
                patient_dictionary[line[0]+'_allele2'] = alleles[1]
            else:
                None


df = pd.DataFrame(all_patient_dictionary).transpose()

df.to_csv('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.alternate.classI.csv')
pickle.dump(all_patient_dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.alternate.classI.p', 'wb'))

# A better dictionary...

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

pickle.dump(dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/Alternate.HLA_classI.p', 'w'))