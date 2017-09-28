import pandas as pd
import scipy.stats as sp
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['legend.numpoints'] = 1
sns.set_style("white")

PATH_TO_GENERATED_FIGURES = '/cellar/users/ramarty/Data/hla_ii/generated_figures/'


def main(category):

    population_frequency(category)

    peptide_class_comparison(category)


# Population statistics
def population_frequency(category):
    patient_affinities = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_affinities.cancer.{0}.csv'.format(category), index_col=0)
    patient_mutations = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.cancer.{0}.csv'.format(category), index_col=0)

    # group mutations by frequency
    counts = list(patient_mutations.sum().unique())
    mutations_by_count = {}
    for count in counts:
        # group all of the mutations with this count in TCGA
        mutations_by_count[count] = list(patient_mutations.sum()[patient_mutations.sum() == count].index)
        # Add column for each patient with the median of their affinities for that mutation count
        patient_affinities[count] = patient_affinities[mutations_by_count[count]].median(axis=1)
    counts.sort()

    grouped = patient_affinities.ix[:, counts].median().reset_index()
    grouped.columns = ['index', 'scores']
    p, rho = sp.spearmanr(grouped.index, grouped.scores)[1], sp.spearmanr(grouped.index, grouped.scores)[0]
    print p, rho

    # would be nice to combine these automatically
    plt.figure(figsize=(16, 2))
    ax = sns.pointplot(x='index', y='scores', data=grouped, order=counts, color="grey", join=False)
    plt.xticks(rotation=45)
    plt.xlabel('Mutations')
    plt.ylabel('Median PHBR')
    plt.title('{0} {1}'.format(p, rho))
    plt.savefig(PATH_TO_GENERATED_FIGURES + 'pop_frequency.scatter.{0}.pdf'.format(category))

    plt.figure(figsize=(20,8))
    sns.heatmap(patient_affinities.ix[:, counts], xticklabels=False, yticklabels=False, vmax=40, cmap=sns.cubehelix_palette(reverse=True, as_cmap=True))
    plt.savefig(PATH_TO_GENERATED_FIGURES + 'pop_frequency.heatmap.{0}.pdf'.format(category))

    plt.figure(figsize=(15.5,2))
    plt.gca().invert_yaxis()
    plt.bar(range(1, len(counts)+1), counts, color='grey')
    plt.xlim(0, 52)
    plt.ylabel('Mutation count in TCGA')
    plt.savefig(PATH_TO_GENERATED_FIGURES + 'pop_frequency.bar.{0}.pdf'.format(category))


# Peptide type
def peptide_class_comparison(category):
    categories = ['oncogenes', 'tsgenes', 'random', 'common', 'viral', 'bacterial']
    df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.cancer.all.csv', index_col=0)
    mutation_counts = pd.DataFrame(df.sum()).reset_index()
    mutation_counts.columns = ['mutation', 'count']
    driver_mutations = list(mutation_counts[mutation_counts['count'] > 10].mutation)
    print len(driver_mutations)
    df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/indels.csv', index_col=0)
    restricted_alleles = [x for x in df.columns if 'DR' in x]

    value_types = []
    for category in categories:
        # restrict the columns to higher frequency mutations
        if category == 'oncogenes' or category == 'tsgenes':
            df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.csv'.format(category), index_col=0)
            df.set_index('mutation', inplace=True)
            app_restricted_space = [x for x in driver_mutations if x in df.index]
            if category == 'all':
                values = get_values_from_df(df.ix[app_restricted_space, :])
            else:
                values = get_values_from_df(df.ix[app_restricted_space, restricted_alleles])
            print category, len(values), len(df.index), len(app_restricted_space)
            value_types.append(values)
        else:
            df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.csv'.format(category), index_col=0)
            df.set_index('mutation', inplace=True)
            if category == 'all':
                values = get_values_from_df(df)
            else:
                values = get_values_from_df(df.ix[:, restricted_alleles])
            value_types.append(values)
            print category, len(values)

    plotting = pd.DataFrame({'category': ['oncogene' for x in value_types[0]] + ['tsgene' for x in value_types[1]] + ['random' for x in value_types[2]] + ['common' for x in value_types[3]] + ['viral' for x in value_types[4]] + ['bacterial' for x in value_types[5]],
                         'best_rank': value_types[0] + value_types[1] + value_types[2] + value_types[3] + value_types[4] + value_types[5]})

    ax = sns.boxplot(x='category', y='best_rank', data=plotting, showfliers=False, color='white')
    ax.grid(False)
    plt.ylabel('Residue Presentation Score')
    plt.xticks(rotation=45)
    plt.xlabel('')
    plt.ylim(0, 100)
    plt.savefig(PATH_TO_GENERATED_FIGURES + 'peptide_class.boxplots.{0}.pdf'.format(category))

    perc_strong, perc_all, total_strong, total_all, total_count = [], [], [], [], []
    for i, category in enumerate(categories):
        total = len(value_types[i])
        greater_than_strong = len(filter(lambda a: a < 2, value_types[i]))
        greater_than_all = len(filter(lambda a: a < 10, value_types[i]))
        total_strong.append(greater_than_strong)
        total_all.append(greater_than_all)
        total_count.append(total)
        print category, float(greater_than_all)/total, float(greater_than_strong)/total
        perc_all.append(float(greater_than_all)/total)
        perc_strong.append(float(greater_than_strong)/total)

    binders = pd.DataFrame({'category': categories, 'Perc_strong': perc_strong, 'Perc_all': perc_all})

    all_binders = pd.DataFrame({'category': list(binders.category) + list(binders.category),
                                'percentage': list(binders.Perc_all) + list(binders.Perc_strong),
                                'binding': ['rank < 10%' for x in list(binders.Perc_all)] + ['rank < 2%' for x in list(binders.Perc_strong)]})

    ax = sns.barplot(x='category', y='percentage', hue='binding', data=all_binders, order=['oncogenes', 'tsgenes', 'random', 'common', 'viral', 'bacterial'], color='white')
    ax.grid(False)
    plt.xticks(rotation=45)
    plt.ylim(0, 0.45)
    plt.ylabel('Fraction of residues with binding peptides')
    plt.xlabel('')
    plt.savefig(PATH_TO_GENERATED_FIGURES + 'peptide_class.percentages.{0}.pdf'.format(category))


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
    main(sys.argv[1])
    sys.exit()