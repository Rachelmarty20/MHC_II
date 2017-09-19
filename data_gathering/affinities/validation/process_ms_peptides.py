import pandas as pd
from Bio import SeqIO
import sys


def main(peptide_length):

    # get proteome reference
    merged = gather_protein_sequences()

    file_names = ['DonorA_DRB10301_DRB11101.csv',  'DonorD_DRB10101_DRB10701.csv',
                  'DonorG_DRB10701_DRB11501_DRB50101.csv', 'DonorB_DRB10401_DRB11301.csv',
                  'DonorE_DRB10101_DRB11101.csv', 'DonorC_DRB10401_DRB11301.csv',
                  'DonorF_DRB10901_DRB11001.csv']

    for donor in ['DonorA', 'DonorB', 'DonorC', 'DonorD', 'DonorE', 'DonorF', 'DonorG']:

        print donor

        file = [x for x in file_names if donor in x][0]

        ms = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/validation/ciudad/raw/{0}'.format(file))
        ms = ms[['Peptide sequence', 'Uniprot AC']]
        ms.columns = ['measured_peptide', 'uniprot_AC']
        ms.measured_peptide = ms.measured_peptide.apply(combine_separated_peptide)
        ms = pd.merge(ms, merged, on='uniprot_AC', how='left')

        def get_combined(x):
            try:
                gene = x[0]
                peptide = x[1]
                aa = peptide[len(peptide)/2]
                searched_position = list(ms[ms.measured_peptide == peptide].sequence)[0].find(peptide)
                residue = searched_position + 1 + (len(peptide) / 2)
                if searched_position == -1:
                    return 'fail'
                else:
                    return '{0}_{1}{2}{3}'.format(gene, aa, residue, aa)
            except:
                return 'fail'
        ms['combined'] = ms[['gene', 'measured_peptide', 'sequence']].apply(get_combined, axis=1)
        processed = ms[ms.combined != 'fail']
        mutations = list(set(processed.combined))

        # output residues
        with open('/cellar/users/ramarty/Data/hla_ii/validation/ciudad/residues/{0}.txt'.format(donor), 'w') as outfile:
            for mutation in mutations:
                outfile.write('{0}\n'.format(mutation))

        processed['peptide_length'] = peptide_length
        processed['sequence_for_affinity'] = processed.loc[:, ['combined', 'sequence', 'peptide_length']].apply(sequence_for_affinity, axis=1)

        with open('/cellar/users/ramarty/Data/hla_ii/validation/ciudad/fasta_files/{0}.{1}.fa'.format(donor, peptide_length), 'w') as outfile:
            for row in processed.iterrows():
                combined = row[1]['combined']
                sequence = row[1]['sequence_for_affinity']
                outfile.write('>{0}\n'.format(combined))
                outfile.write('{0}\n'.format(sequence))

    # random
    merged = gather_protein_sequences_all()

    for donor in ['DonorA', 'DonorB', 'DonorC', 'DonorD', 'DonorE', 'DonorF', 'DonorG']:

        print donor

        mutations = [x.strip() for x in open('/cellar/users/ramarty/Data/hla_ii/presentation/residues/random.txt').readlines()][:1000]
        with open('/cellar/users/ramarty/Data/hla_ii/validation/ciudad/residues/{0}.random.txt'.format(donor), 'w') as outfile:
            for mutation in mutations:
                outfile.write('{0}\n'.format(mutation))

        peptides, mutations_used = generate_peptides(mutations, merged, peptide_length)

        with open('/cellar/users/ramarty/Data/hla_ii/validation/ciudad/fasta_files/{0}.random.{1}.fa'.format(donor, peptide_length), 'w') as outfile:
            for mutation, sequence in zip(mutations_used, peptides):
                outfile.write('>{0}\n'.format(mutation))
                outfile.write('{0}\n'.format(sequence))




def gather_protein_sequences():
    # Ensemble (with sequences)
    proteins, genes, sequences = [], [], []
    fasta_sequences = SeqIO.parse(open("/cellar/users/ramarty/Data/hla_ii/references/Homo_sapiens.GRCh38.pep.all.fa"),'fasta')
    for fasta in fasta_sequences:
        gene, protein, sequence = fasta.description.split('gene_symbol:')[1].split(' ')[0], fasta.id, fasta.seq.tostring()
        proteins.append(protein)
        genes.append(gene)
        sequences.append(sequence)
    ensemble = pd.DataFrame({'protein': proteins, 'gene': genes, 'sequence': sequences})
    ensemble['protein'] = ensemble['protein'].map(lambda x: str(x)[:-2])

    # Uniprot (mapping)
    uniprot_import = pd.read_csv('/cellar/users/ramarty/Data/hla_ii/references/HUMAN_9606_idmapping_selected.tab', sep='\t', header=None)
    proteins, uniprot_ACs = [], []
    for i in range(len(uniprot_import)):
        try:
            proteins.append(uniprot_import.ix[i, 20].split(';')[0])
            uniprot_ACs.append(uniprot_import.ix[i, 0])
        except:
            None
    uniprot = pd.DataFrame({'protein': proteins, 'uniprot_AC': uniprot_ACs})
    # Drop all rows without a Uniprot AC
    merged = pd.merge(ensemble, uniprot, on='protein', how='left')
    merged = merged.drop_duplicates()
    merged = merged.dropna()
    # Merge to make single, connonical entry per gene
    merged['length'] = merged['sequence'].apply(get_length)
    merged = merged.sort_values('length', ascending=False)
    merged = merged.drop_duplicates(subset='protein')
    return merged


def gather_protein_sequences_all():
    # Ensemble (with sequences)
    proteins, genes, sequences = [], [], []
    fasta_sequences = SeqIO.parse(open("/cellar/users/ramarty/Data/hla_ii/references/Homo_sapiens.GRCh38.pep.all.fa"),'fasta')
    for fasta in fasta_sequences:
        gene, protein, sequence = fasta.description.split('gene_symbol:')[1].split(' ')[0], fasta.id, fasta.seq.tostring()
        proteins.append(protein)
        genes.append(gene)
        sequences.append(sequence)
    ensemble = pd.DataFrame({'protein': proteins, 'gene': genes, 'sequence': sequences})
    ensemble['protein'] = ensemble['protein'].map(lambda x: str(x)[:-2])
    # Add missing gene
    ensemble.loc[len(ensemble)+1] = ['CRIPAK', '_', 'MHEPSLCANVECPPAHTCPCGVPACSCAHVECPPAHTCRCGVPACSHMPMWSARLLTRAHVECPPAHTRVHVECPPAHVPMWSAHLLTCADVECHLLTHVPMWSARLLTCPCGVPACSHVPMRSARLLTRAHAECPPAHTCPCGVPACSHVPMRSARLLTRADVECPPAHTCPCGVPACSHVPTWSARLITRAHVECSPAHTCRCGVPACSHVPMWSVRLLTRADAECPPAHTCRCGVPACSHVPMWSARLLTCRCGVPACSHVPMWSARLLTCRCGVPACSHVPMWSARLLTRAHVECPPAHTCRRGVPACSRAHMECPPAHTCHCGVPACSHTCRCGVPACSHVPMWSARLLTRAHVECPPAHTRAHVECPPAHTCPCGVPACSHTCPCGVPACSHKALAWWFCRFPVLPAESDAVTVHSTHGGFLIRFYVKDPFYISLHLEIT']
    return ensemble


def get_length(x):
    return len(x)


def combine_separated_peptide(x):
    if len(x.split(' ')) > 1:
        return ''.join(x.split(' '))
    else:
        return x


def sequence_for_affinity(x):
    combined = x[0]
    sequence = x[1]
    peptide_length = x[2]
    residue = combined.split('_')[1]
    position = int(residue[1:len(residue)-1]) - 1
    old_aa = residue[0]
    new_aa = residue[-1:]
    if len(sequence) >= position and sequence[position] == old_aa:
        mutated_sequence = sequence[:position] + new_aa + sequence[position+1:]
        if position > (peptide_length-2):
            seq_for_affinity = mutated_sequence[position-(peptide_length-1):position+(peptide_length)]
        else:
            seq_for_affinity = mutated_sequence[:position+peptide_length]
        return seq_for_affinity
    else:
        return 'fail'


def generate_peptides(mutations, merged, peptide_length):
    peptides, mutations_used = [], []
    for mutation in mutations:
        gene = mutation.split('_')[0]
        sequences = list(merged[merged.gene == gene].sequence)
        for i, sequence in enumerate(sequences):
            residue = mutation.split('_')[1]
            position = int(residue[1:len(residue)-1]) - 1
            old_aa = residue[0]
            new_aa = residue[-1:]
            if len(sequence) > position and sequence[position] == old_aa:
                mutated_sequence = sequence[:position] + new_aa + sequence[position+1:]
                if position > (peptide_length - 2):
                    seq_for_affinity = mutated_sequence[position-(peptide_length-1):position+(peptide_length)]
                else:
                    seq_for_affinity = mutated_sequence[:position+(peptide_length)]
                peptides.append(seq_for_affinity)
                mutations_used.append(mutation)
                break
            else:
                if i+1 == len(sequences):
                    #print mutation
                    continue
    return peptides, mutations_used

###########################################  Main Method  #####################################

if __name__ == "__main__":
    if len(sys.argv) != 1:
        sys.exit()
    main(int(sys.argv[1]))
    sys.exit()