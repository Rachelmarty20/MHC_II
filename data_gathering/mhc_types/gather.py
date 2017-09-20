import pandas as pd
import cPickle as pickle
import samples

# get TCGA types
barcodes = samples.get_barcodes()

all_patient_dictionary = {}
for i, barcode in enumerate(barcodes[:2000]):

    patient_dictionary = {}
    directory = '/nrnb/users/ramarty/TCGA/exomes/{0}/hlaHD/sampleID/result/'.format(barcode)

    try:
        lines = [x.split('\t') for x in open('{0}sampleID_final.result.txt'.format(directory)).readlines()[:8]]

        patient_dictionary = {}
        for line in lines:
            if line[0] in ['DRB1']:
                if line[2].strip() == '-':
                    x = line[1]
                    alleles = [line[0]+'_'+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1]]*2
                else:
                    alleles = [line[0]+'_'+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1] for x in line[1:3]]
                patient_dictionary[line[0]+'_allele1'] = alleles[0]
                patient_dictionary[line[0]+'_allele2'] = alleles[1]
            # combinations
            elif line[0] in ['DQA1', 'DQB1', 'DPA1', 'DPB1']:
                if line[2].strip() == '-':
                    x = line[1]
                    alleles = [line[0]+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1]]*2
                else:
                    alleles = [line[0]+x.split('*')[1].split(':')[0]+x.split('*')[1].split(':')[1] for x in line[1:3]]
                patient_dictionary[line[0]+'_allele1'] = alleles[0]
                patient_dictionary[line[0]+'_allele2'] = alleles[1]
                all_patient_dictionary[barcode] = patient_dictionary
            else:
                None

    except:
        print 'skip' + barcode

df = pd.DataFrame(all_patient_dictionary).transpose()

df.to_csv('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.tcga.csv')
pickle.dump(all_patient_dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.tcga.p', 'wb'))


# HLA-DPA10103-DPB10101
# DRB1_1302