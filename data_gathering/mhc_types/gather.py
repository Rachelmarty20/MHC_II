import pandas as pd
import cPickle as pickle
import samples


# genes to collect
genes = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'DRB3', 'DRB4']

# get TCGA types
barcodes = samples.get_barcodes()

all_patient_dictionary = {}
for i, barcode in enumerate(barcodes[:5]):
    print barcode

    patient_dictionary = {}
    directory = '/nrnb/users/ramarty/TCGA/exomes/{0}/hlaHD/sampleID/result/'.format(barcode)

    for gene in genes:
        print directory + 'sampleID_{0}.est.txt'.format(gene)
        patient_dictionary[gene] = [':'.join(x.split(':')[:2]) for x in open(directory + 'sampleID_{0}.est.txt'.format(gene)).readlines()[2].split('\t')[:2]]

    all_patient_dictionary[barcode] = patient_dictionary

df = pd.DataFrame(all_patient_dictionary).transpose()

df.to_csv('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.tcga.csv')
pickle.dump(df, open('/cellar/users/ramarty/Data/hla_ii/hla_types/hla_types.tcga.p', 'wb'))