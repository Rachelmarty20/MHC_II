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
PATH_TO_DATA = '/cellar/users/ramarty/Data/hla/git_data/'

def main(category):

    patient_affinities = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_affinities.cancer.{0}.csv'.format(category), index_col=0)
    patient_mutations = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.cancer.{0}.csv'.format(category), index_col=0)

    population_frequency(category, patient_affinities, patient_mutations)

    peptide_class_comparison(category, patient_affinities, patient_mutations)

    heatmap_overview(category, patient_affinities, patient_mutations)

# Heatmap
def heatmap_overview(category, patient_affinities, patient_mutations):

    all_mutations = list(patient_affinities.columns)
    oncogenes = [x.strip() for x in open(PATH_TO_DATA + 'data/onco_genes.txt').readlines()]
    tsgenes = [x.strip() for x in open(PATH_TO_DATA + 'data/tumor_suppressor_genes.txt').readlines()]
    onco_muts = [x for x in all_mutations if x.split('_')[0] in oncogenes]
    ts_muts = [x for x in all_mutations if x.split('_')[0] in tsgenes]
    missense_muts = [x for x in all_mutations if len(x.split('_')) == 2]
    indel_muts = [x for x in all_mutations if len(x.split('_')) == 3]

    clinical = pd.read_csv(PATH_TO_DATA + 'data/clinical/ancestory.csv', index_col=0)
    asian = pd.Series([x for x in list(clinical[clinical.race.isin(['ASIAN'])].index) if x in list(patient_affinities.index)]).sample(100)
    black = pd.Series([x for x in list(clinical[clinical.race.isin(['BLACK OR AFRICAN AMERICAN'])].index) if x in list(patient_affinities.index)]).sample(100)
    white = pd.Series([x for x in list(clinical[clinical.race.isin(['WHITE'])].index) if x in list(patient_affinities.index)]).sample(100)

    restricted_patients = list(asian) + list(black) + list(white)
    patient_affinities_small = patient_affinities.ix[restricted_patients, :]
    patient_mutations_small = patient_mutations.ix[restricted_patients, :]

    # add mutation type
    mutation_df = pd.DataFrame({'Classification': ['Missense' for x in missense_muts if x in all_mutations] + ['Indel' for x in indel_muts if x in all_mutations],
                                     'Mutation': [x for x in missense_muts if x in all_mutations] + [x for x in indel_muts if x in all_mutations]})
    cmap = sns.color_palette("Paired", 2)
    def map_colors(x):
        if x == 'Missense':
            return cmap[0]
        elif x == 'Indel':
            return cmap[1]
    mutation_df['Color'] = mutation_df['Classification'].apply(map_colors)
    # add oncogene or tumor suppressor
    gene_df = pd.DataFrame({'Classification': ['Onco' for x in onco_muts if x in all_mutations] + ['TS' for x in ts_muts if x in all_mutations],
                                     'Mutation': [x for x in onco_muts if x in all_mutations] + [x for x in ts_muts if x in all_mutations]})
    cmap = sns.color_palette("Paired", 4)
    def map_colors(x):
        if x == 'Onco':
            return cmap[2]
        elif x == 'TS':
            return cmap[3]
    gene_df['Color'] = gene_df['Classification'].apply(map_colors)
    sns.palplot(sns.color_palette("Paired", 4)[2:])

    # add patient tissue - top 7 colors and then an 'other' color
    patient_df = pd.read_csv(PATH_TO_DATA + 'data/patient_tissues.csv', index_col=0)
    patient_df = patient_df.ix[restricted_patients]
    large_tissues = list(patient_df.Tissue.value_counts().index[:7])
    cmap = sns.color_palette("Set2", 10)
    def map_colors(x):
        if x == large_tissues[0]:
            return cmap[0]
        elif x == large_tissues[1]:
            return cmap[1]
        elif x == large_tissues[2]:
            return cmap[2]
        elif x == large_tissues[3]:
            return cmap[3]
        elif x == large_tissues[4]:
            return cmap[4]
        elif x == large_tissues[5]:
            return cmap[5]
        elif x == large_tissues[6]:
            return cmap[6]
        else:
            # maybe change 'other' to white?
            return (1, 1, 1)
    patient_df['Color'] = patient_df['Tissue'].apply(map_colors)
    print large_tissues
    sns.palplot(sns.color_palette("Set2", 8))

    # add ethnicity to patient_df
    clinical = pd.read_csv(PATH_TO_DATA + 'data/clinical/ancestory.csv', index_col=0)
    patient_df = pd.merge(patient_df, clinical[['race']], left_index=True, right_index=True, how='left')
    patient_df['race'] = list(patient_df['race'].fillna('OTHER'))
    cmap = sns.color_palette("BrBG", 4)
    def map_colors(x):
        if x == 'WHITE':
            return cmap[1]
        elif x == 'BLACK OR AFRICAN AMERICAN':
            return cmap[0]
        elif x == 'ASIAN':
            return cmap[3]
        else:
            return cmap[2]
    patient_df['Ethnicity_Color'] = patient_df['race'].apply(map_colors)
    sns.palplot([sns.color_palette("BrBG", 4)[1],sns.color_palette("BrBG", 4)[0],sns.color_palette("BrBG", 4)[3]])

    # add immune infiltration
    hot_patients = [x.strip() for x in open(PATH_TO_DATA + 'data/clinical/immune_infiltration/hot_patients.txt').readlines() if x.strip() in restricted_patients]
    cold_patients = [x.strip() for x in open(PATH_TO_DATA + 'data/clinical/immune_infiltration/cold_patients.txt').readlines() if x.strip() in restricted_patients]
    warm_patients = [x.strip() for x in open(PATH_TO_DATA + 'data/clinical/immune_infiltration/warm_patients.txt').readlines() if x.strip() in restricted_patients]
    cool_patients = [x.strip() for x in open(PATH_TO_DATA + 'data/clinical/immune_infiltration/cool_patients.txt').readlines() if x.strip() in restricted_patients]
    immune_infil = []
    for patient in list(patient_df.index):
        if patient in hot_patients:
            immune_infil.append('hot')
        elif patient in cold_patients:
            immune_infil.append('cold')
        elif patient in cool_patients:
            immune_infil.append('cool')
        elif patient in warm_patients:
            immune_infil.append('warm')
        else:
            immune_infil.append('neither')
    patient_df['Immune_infil'] = immune_infil
    cmap = sns.color_palette("RdBu", 8)
    def map_colors(x):
        if x == 'hot':
            return cmap[0]
        elif x == 'warm':
            return cmap[1]
        elif x == 'neither':
            return (1,1,1)
        elif x == 'cool':
            return cmap[6]
        elif x == 'cold':
            return cmap[7]
    patient_df['Immune_Color'] = patient_df['Immune_infil'].apply(map_colors)
    sns.palplot(sns.color_palette("RdBu", 8)[:2] + sns.color_palette("RdBu", 8)[-2:])

    plt.figure(figsize=(12,12))
    sns.clustermap(patient_affinities_small.ix[list(patient_df.index), list(mutation_df.Mutation)],#, xticklabels=False, yticklabels=False,
               row_colors=[list(patient_df.Ethnicity_Color), list(patient_df.Immune_Color), list(patient_df.Color)],
               col_colors=[list(mutation_df.Color), list(gene_df.Color)], vmax=40, cmap=sns.cubehelix_palette(reverse=True, as_cmap=True))
    plt.savefig(PATH_TO_GENERATED_FIGURES + 'clustermap.PHBR.pdf')


# Population statistics
def population_frequency(category, patient_affinities, patient_mutations):
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

    # would be nice to combine these automatically
    plt.figure(figsize=(16, 2))
    ax = sns.pointplot(x='index', y='scores', data=grouped, order=counts, color="grey", join=False)
    plt.xticks(rotation=45)
    plt.xlabel('Mutations')
    plt.ylabel('Median PHBR')
    plt.title('{0} {1}'.format(p, rho))
    plt.savefig(PATH_TO_GENERATED_FIGURES + '{0}/pop_frequency.scatter.pdf'.format(category))
    plt.clf()

    plt.figure(figsize=(20,8))
    sns.heatmap(patient_affinities.ix[:, counts], xticklabels=False, yticklabels=False, vmax=40, cmap=sns.cubehelix_palette(reverse=True, as_cmap=True))
    plt.savefig(PATH_TO_GENERATED_FIGURES + '{0}/pop_frequency.heatmap.pdf'.format(category))
    plt.clf()

    plt.figure(figsize=(15.5,2))
    plt.gca().invert_yaxis()
    plt.bar(range(1, len(counts)+1), counts, color='grey')
    plt.xlim(0, 52)
    plt.ylabel('Mutation count in TCGA')
    plt.savefig(PATH_TO_GENERATED_FIGURES + '{0}/pop_frequency.bar.pdf'.format(category))
    plt.clf()

# Peptide type
def peptide_class_comparison(category, patient_affinities, patient_mutations):
    categories = ['oncogenes', 'tsgenes', 'random', 'common', 'viral', 'bacterial']
    df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/clean_matrices/patient_mutations.cancer.all.csv', index_col=0)
    mutation_counts = pd.DataFrame(df.sum()).reset_index()
    mutation_counts.columns = ['mutation', 'count']
    driver_mutations = list(mutation_counts[mutation_counts['count'] > 10].mutation)
    df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/indels.csv', index_col=0)
    restricted_alleles = [x for x in df.columns if 'DR' in x]

    value_types = []
    for cat in categories:
        # restrict the columns to higher frequency mutations
        if cat == 'oncogenes' or cat == 'tsgenes':
            df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.csv'.format(cat), index_col=0)
            df.set_index('mutation', inplace=True)
            app_restricted_space = [x for x in driver_mutations if x in df.index]
            if cat == 'all':
                values = get_values_from_df(df.ix[app_restricted_space, :])
            else:
                values = get_values_from_df(df.ix[app_restricted_space, restricted_alleles])
            value_types.append(values)
        else:
            df = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/presentation/allele_matrices/{0}.csv'.format(cat), index_col=0)
            df.set_index('mutation', inplace=True)
            if cat == 'all':
                values = get_values_from_df(df)
            else:
                values = get_values_from_df(df.ix[:, restricted_alleles])
            value_types.append(values)

    plotting = pd.DataFrame({'kind': ['oncogene' for x in value_types[0]] + ['tsgene' for x in value_types[1]] + ['random' for x in value_types[2]] + ['common' for x in value_types[3]] + ['viral' for x in value_types[4]] + ['bacterial' for x in value_types[5]],
                         'best_rank': value_types[0] + value_types[1] + value_types[2] + value_types[3] + value_types[4] + value_types[5]})

    plt.figure(figsize=(8,6))
    ax = sns.boxplot(x='kind', y='best_rank', data=plotting, showfliers=False, color='white')
    ax.grid(False)
    plt.ylabel('Residue Presentation Score')
    plt.xticks(rotation=45)
    plt.xlabel('')
    plt.ylim(0, 100)
    plt.savefig(PATH_TO_GENERATED_FIGURES + '{0}/peptide_class.boxplots.pdf'.format(category))
    plt.clf()

    perc_strong, perc_all, total_strong, total_all, total_count = [], [], [], [], []
    for i, cat in enumerate(categories):
        total = len(value_types[i])
        greater_than_strong = len(filter(lambda a: a < 2, value_types[i]))
        greater_than_all = len(filter(lambda a: a < 10, value_types[i]))
        total_strong.append(greater_than_strong)
        total_all.append(greater_than_all)
        total_count.append(total)
        perc_all.append(float(greater_than_all)/total)
        perc_strong.append(float(greater_than_strong)/total)

    binders = pd.DataFrame({'kind': categories, 'Perc_strong': perc_strong, 'Perc_all': perc_all})

    all_binders = pd.DataFrame({'kind': list(binders.kind) + list(binders.kind),
                                'percentage': list(binders.Perc_all) + list(binders.Perc_strong),
                                'binding': ['rank < 10%' for x in list(binders.Perc_all)] + ['rank < 2%' for x in list(binders.Perc_strong)]})

    plt.figure(figsize=(8,6))
    ax = sns.barplot(x='kind', y='percentage', hue='binding', data=all_binders, order=['oncogenes', 'tsgenes', 'random', 'common', 'viral', 'bacterial'], color='white')
    ax.grid(False)
    plt.xticks(rotation=45)
    plt.ylim(0, 0.45)
    plt.ylabel('Fraction of residues with binding peptides')
    plt.xlabel('')
    plt.savefig(PATH_TO_GENERATED_FIGURES + '{0}/peptide_class.percentages.pdf'.format(category))
    plt.clf()

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
