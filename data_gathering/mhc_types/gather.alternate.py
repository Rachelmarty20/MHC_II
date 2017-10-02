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
            if line[0] in ['DRB1']:
                if line[1] == 'Not typed':
                    alleles = ['-', '-']
                elif line[2].strip() == '-':
                    x = line[1]
                    alleles = [line[0]+'_'+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1]]*2
                else:
                    alleles = [line[0]+'_'+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1] for x in line[1:3]]
                patient_dictionary[line[0]+'_allele1'] = alleles[0]
                patient_dictionary[line[0]+'_allele2'] = alleles[1]
            # combinations
            elif line[0] in ['DQA1', 'DQB1', 'DPA1', 'DPB1']:
                if line[1] == 'Not typed':
                    alleles = ['-', '-']
                elif line[2].strip() == '-':
                    x = line[1]
                    alleles = [line[0]+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1]]*2
                else:
                    alleles = [line[0]+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1] for x in line[1:3]]
                patient_dictionary[line[0]+'_allele1'] = alleles[0]
                patient_dictionary[line[0]+'_allele2'] = alleles[1]
                all_patient_dictionary[sample] = patient_dictionary
            else:
                None


df = pd.DataFrame(all_patient_dictionary).transpose()

df.to_csv('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.alternate.csv')
pickle.dump(all_patient_dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.alternate.p', 'wb'))

# A better dictionary...

patient_types_df = df
patients = list(patient_types_df.index)

dictionary = {}
for patient in patients:
    patient_dictionary = {}
    try:
        # DR - adding two of each for consistency sake with other genes
        alleles = list(patient_types_df.ix[patient][['DRB1_allele1', 'DRB1_allele2']])
        alleles.extend(list(patient_types_df.ix[patient][['DRB1_allele1', 'DRB1_allele2']]))
        patient_dictionary['DR'] = alleles

        # DP
        alleles = []
        DA = list(patient_types_df.ix[patient][['DPA1_allele1', 'DPA1_allele2']])
        DB = list(patient_types_df.ix[patient][['DPB1_allele1', 'DPB1_allele2']])
        if ('-' not in DA) & ('-' not in DB):
            for a in DA:
                for b in DB:
                    alleles.append('HLA-{0}-{1}'.format(a.strip(), b.strip()))
        patient_dictionary['DP'] = alleles

        # DQ
        alleles = []
        DA = list(patient_types_df.ix[patient][['DQA1_allele1', 'DQA1_allele2']])
        DB = list(patient_types_df.ix[patient][['DQB1_allele1', 'DQB1_allele2']])
        if ('-' not in DA) & ('-' not in DB):
            for a in DA:
                for b in DB:
                    alleles.append('HLA-{0}-{1}'.format(a.strip(), b.strip()))
        patient_dictionary['DQ'] = alleles

        dictionary[patient] = patient_dictionary
    except:
        # weeds out the patients with the failed typing
        print patient

pickle.dump(dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/Alternate.HLA_classII.p', 'w'))