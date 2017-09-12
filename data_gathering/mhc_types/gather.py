import pandas as pd
import cPickle as pickle
import samples

# genes to collect
genes = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'DRB3', 'DRB4']

# get TCGA types
barcodes = samples.get_barcodes()

all_patient_dictionary = {}
for i, barcode in enumerate(barcodes[:2000]):

    patient_dictionary = {}
    directory = '/nrnb/users/ramarty/TCGA/exomes/{0}/hlaHD/sampleID/result/'.format(barcode)

    try:
        open(directory + 'sampleID_A.est.txt')

        for gene in genes:

            try:
                patient_dictionary[gene] = [':'.join(x.split(':')[:2]) for x in open(directory + 'sampleID_{0}.est.txt'.format(gene)).readlines()[2].split('\t')[:2]]
            except:
                patient_dictionary[gene] = ['-', '-']

        all_patient_dictionary[barcode] = patient_dictionary

    except:
        print 'skip' + barcode

df = pd.DataFrame(all_patient_dictionary).transpose()

df.to_csv('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.tcga.csv')
pickle.dump(all_patient_dictionary, open('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.tcga.p', 'wb'))