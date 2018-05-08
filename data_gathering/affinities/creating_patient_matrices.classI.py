import pandas as pd
import cPickle as pickle
import sys

# Update to include options for


def main(category, population, subset, condition):

    patient_dictionary = pickle.load(open('/cellar/users/ramarty/Data/hla_ii/hla_types/clean/{0}.ClassI.{1}.p'.format(population,
                                                                                                                       subset)))

    df = pd.read_csv('/cellar/users/ramarty/Data/hla/residue/matrices/{0}/class_i/{1}.all.max.txt'.format(condition,
                                                                                               category), index_col=0)
    patients_used = []
    for patient in patient_dictionary.keys():

        df[patient] = df[patient_dictionary[patient]].apply(PHBR, axis=1)
        patients_used.append(patient)

    df[patients_used].to_csv('/cellar/users/ramarty/Data/hla_ii/presentation/patient_matrices/{0}.{1}.{2}.{3}.ClassI.csv'.format(category,
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
