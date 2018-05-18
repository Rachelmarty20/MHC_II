import pandas as pd
import cPickle as pickle
import sys


def main(iterations):

    df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.cancer.TCGA.inclusive.mut.csv',
                 index_col=0)
    mutation_counts = pd.DataFrame(df.sum()).reset_index()
    mutation_counts.columns = ['mutation', 'count']
    driver_mutations = list(mutation_counts[mutation_counts['count'] > 10].mutation)
    #driver_mutations = list(mutation_counts.mutation)
    print len(driver_mutations)
    categories = ['oncogenes', 'tsgenes', 'random0', 'germline', 'viral', 'bacterial']
    value_types = []
    for category in categories:
        # restrict the columns to higher frequency mutations
        if category == 'oncogenes' or category == 'tsgenes':
            df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/patient_matrices/{0}.TCGA.inclusive.mut.csv'.format(category),
                             index_col=0)
            #df.set_index('mutation', inplace=True)
            app_restricted_space = [x for x in driver_mutations if x in df.index]
            values = get_values_from_df(df.ix[app_restricted_space, :])
            print category, len(values), len(df.index), len(app_restricted_space)
            value_types.append(values)
        else:
            try:
                df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/patient_matrices/{0}.TCGA.inclusive.mut.csv'.format(category),
                                 index_col=0)
            except:
                df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/patient_matrices/{0}.all.csv'.format(category),
                                 index_col=0)
            #df.set_index('mutation', inplace=True)
            values = get_values_from_df(df)
            value_types.append(values)
            print category, len(values)
    categories = ['Oncogenes', 'TSgenes', 'Random', 'Germline', 'Viral', 'Bacterial']
    plotting = pd.DataFrame({'category': ['Oncogenes' for x in value_types[0]] + ['TSgenes' for x in value_types[1]] + ['Random' for x in value_types[2]] + ['Germline' for x in value_types[3]] + ['Viral' for x in value_types[4]] + ['Bacterial' for x in value_types[5]],
                             'PHBR': value_types[0] + value_types[1] + value_types[2] + value_types[3] + value_types[4] + value_types[5]})

    fraction = {}
    for cat in categories:
        fraction[cat] = []
    for cat in categories:
        print cat
        tmp = plotting[plotting.category == cat].PHBR
        for i in range(iterations):
            single_trajectory = []
            tmp_boot = tmp.sample(len(tmp), replace = True)
            for x in np.arange(0, 31, 1):
                single_trajectory.append(mean(tmp_boot < x))
            fraction[cat].append(single_trajectory)

    # calculate confidence intervals at a cutoff of 6
    CI_lb, CI_ub, median = [], [], []
    for cat in categories:
        distribution = []
        for i in range(len(fraction[cat])):
            distribution.append(fraction[cat][i][6])
        CI_lb.append(pd.Series(distribution).quantile(0.025))
        median.append(pd.Series(distribution).quantile(0.5))
        CI_ub.append(pd.Series(distribution).quantile(0.975))
    results = pd.DataFrame({'Category': categories,
                              'Median': median,
                              'Lower_bound': CI_lb,
                              'Upper_bound': CI_ub})
    results.to_csv('/cellar/users/ramarty/Data/hla_ii/generated_data/peptide_comparison/confidence_intervals.threshold_6.iterations_{0}.csv'.format(iterations))



def get_values_from_df(df):
    values = []
    for a in df.values:
        values.extend(a)
    return values


###########################################  Main Method  #####################################

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print 'Wrong number of arguments.'
        sys.exit()
    main(int(sys.argv[1]))
    sys.exit()
