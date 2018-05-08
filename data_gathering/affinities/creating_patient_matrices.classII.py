import pandas as pd
import cPickle as pickle
import sys

# Update to include options for


def main(category, population, subset, condition):

    patient_dictionary = pickle.load(open('/cellar/users/ramarty/Data/hla_ii/hla_types/clean/{0}.ClassII.{1}.p'.format(population,
                                                                                                                       subset)))

    if condition == 'mut':
        df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.csv'.format(category), index_col=0)
    else:
        df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.wt.csv'.format(category), index_col=0)
    patients_used = []
    for patient in patient_dictionary.keys():
        patient_alleles = []
        for gene in ['DR', 'DP', 'DQ']:
            patient_alleles.extend(patient_dictionary[patient][gene])
        try:
            if len(patient_alleles) == 12:
                df[patient] = df[patient_alleles].apply(PHBR, axis=1)
                patients_used.append(patient)
        except:
            print patient
    df.index = df['mutation']
    df[patients_used].to_csv('/cellar/users/ramarty/Data/hla_ii/presentation/patient_matrices/{0}.{1}.{2}.{3}.csv'.format(category,
                                                                                                                              population,
                                                                                                                              subset,
                                                                                                                              condition))

def PHBR(x):
    number_of_alleles = len(x)
    s = 0
    for i in range(number_of_alleles):
        s += 1/float(x[i])
    return number_of_alleles / float(s)


###########################################  Main Method  #####################################

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print 'Wrong number of arguments.'
        sys.exit()
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    sys.exit()
