import csv
import math
import os

import numpy as np
import pandas as pd
import config as cfg
import json
import matplotlib.pyplot as plt
import seaborn as sns

terms_names = {

    '00078': 'pre-molten globule',
    '00077': 'molten globule',
    '00051': 'disorder to molten globule',
    '00052': 'disorder to pre-molten globule',
    '00053': 'molten globule to order',
    '00055': 'pre-molten globule to order',
    '00058': 'Order to pre-molten globule',
    '00060': 'Pre-molten globule to molten globule',
    '00063': 'Protein binding',
    '00065': 'DNA binding',
    '00066': 'RNA binding',
    '00071': 'Lipid binding',
    '00073': 'Metal ion binding',
    '00074': 'Small molecule binding',
    '00002': 'Flexible linker/spacer',
    '00003': 'Entropic bristle',
    '00005': 'Entropic spring',
    '00006': 'Structural mortar',
    '00007': 'Self-transport through channel',
    '00009': 'Assembler',
    '00011': 'Targeting',
    '00012': 'Tethering',
    '00014': 'Neutralization of toxic molecules',
    '00015': 'Metal binding',
    '00016': 'Water storage',
    '00019': 'Cis-regulatory element',
    '00020': 'Disassembler',
    '00021': 'Activator',
    '00022': 'DNA bending',
    '00023': 'DNA unwinding',
    '00025': 'Regulation of phosphorylation',
    '00026': 'Regulation of acetylation',
    '00027': 'Regulation of methylation',
    '00028': 'Regulation of glycosylation',
    '00029': 'Regulation of ubiquitination',
    '00032': 'Regulation of palmitoylation',
    '00033': 'Regulation of limited proteolysis',
    '00034': 'Regulation of ADP-ribosylation',
    '00037': 'Space filling',
    '00038': 'Entropic exclusion',
    '00041': 'Liquid-liquid phase separation',
    '00042': 'Hydrogel',
    '00043': 'Aggregate',
    '00046': 'Amyloid',
    '00047': 'Prion',

}


def create_protein_terms_features_table(terms, features):
    # for term in terms:
        # parse_file_protein_terms(term)

    for feature in features:
            # add features to files
            parse_file_proteins_features(terms, feature)
        #     generate_heatmap_dataset(term, feature)
        # generate_heatmaps(feature)
        # I could erase all generated files with a command here at the end


def parse_file_protein_terms(term):

    with open(cfg.data['data'] + 'disprot_regions.json') as json_file:
        data = json.load(json_file)
        dict_terms = {}
        df = pd.DataFrame()
        the_term_is_in_the_whole_file = False

        for el in data['data']:
            the_term_is_in_the_current_row = False

            sequence_length = len(el['sequence'])

            for region in el['regions']:
                if region['term_id'] == term:
                    the_term_is_in_the_whole_file = True
                    the_term_is_in_the_current_row = True

                    for i in range(region['start'] - 1, region['end'] - 1):
                        if el['acc'] in dict_terms:
                            dict_terms[el['acc']][i] = 1
                        else:
                            dict_terms[el['acc']] = [0] * sequence_length
            if the_term_is_in_the_current_row:
                df_el = pd.DataFrame({
                    'residue' : np.arange(1, sequence_length + 1, 1).tolist(),
                    'acc': el['acc'],
                    'term_code': terms_names[term],
                    'term': dict_terms[el['acc']]
                })
                df = df.append(df_el)

        if the_term_is_in_the_whole_file == False:
            print(term + ' not found')

        file = cfg.data['data'] + 'proteins_terms_features/terms_features.tsv'
        df.to_csv(file, mode='a',  header=(not os.path.exists(file)), sep='\t', index=False)




def parse_file_proteins_features(terms, feature):
    full_feature_column = []

    for term in terms:

        df_features = pd.read_csv(cfg.data['data'] + 'mobidb_features.tsv', sep='\t')
        print(feature)

        # calculate features for each acc
        terms_features_file = cfg.data['data'] + 'proteins_terms_features/terms_features.tsv'
        df_terms_features = pd.read_csv(terms_features_file, sep= '\t')
        accessions = df_terms_features.acc.unique()
        for acc in accessions:

            sequence_length = df_terms_features.loc[(df_terms_features.acc == acc) & (df_terms_features.term_code == terms_names[term]) ]
            sequence_length = len(sequence_length.index)
            if sequence_length == 0:
                continue
            print(acc)
            print(sequence_length)

            regions = df_features.loc[(df_features.acc == acc) & (df_features.feature == feature)]['start..end'].values
            if len(regions) > 0:
                regions = regions[0].split(',')
                print(regions)

            # write features to file
            # will be zero if feature is missing or if acc is missing in the mobidb file
            feature_column = [0] * sequence_length
            for region in regions:
                start = region.split('..')[0]
                end = region.split('..')[1]
                for i in range(int(start) - 1, int(end) - 1):
                    print(i)
                    feature_column[i] = 1
            full_feature_column = full_feature_column + feature_column

    df_terms_features[feature] = full_feature_column
    df_terms_features.to_csv(terms_features_file, sep='\t')


def calculate_jaccard_index(term, feature):
    dict_jaccard = {term + '-' + feature: []}
    term_dir = cfg.data['data'] + 'proteins_terms_features/' + term + '/'

    for filename in os.listdir(term_dir):
        df_file = pd.read_csv(term_dir + filename, sep='\t', index_col=[0])
        intersection = df_file.loc[(df_file.term == 1) & (df_file[feature] == 1)]
        union = df_file.loc[(df_file.term == 1) | (df_file[feature] == 1)]
        jaccard = 0
        if len(intersection) != 0:
            jaccard = round(len(intersection) / len(union), 2)
            # print(jaccard)
        dict_jaccard[term + '-' + feature].append(round(jaccard, 1))


    row = []
    # how many times each step happens
    total_occurrences = []

    for j in np.arange(0, 1, 0.1):
        jaccard_occurrence = dict_jaccard[term + '-' + feature].count(round(j, 1))
        total_occurrences.append(jaccard_occurrence)

    print(total_occurrences)
    for x in range(10):
        row.append(term)  # TODO: add quantity of proteins for each term
        row.append(x + 1)
        # normalizing the occurrences to the range 0 1
        row.append(total_occurrences[x] / sum(total_occurrences))

    row_splitted = [row[i:i + 3] for i in range(0, len(row), 3)]

    return row_splitted


def generate_heatmap_dataset(term, feature):
    rows = calculate_jaccard_index(term, feature)
    headers = ['term', 'overlap', 'proteins', 'feature']
    file_exists = os.path.isfile(cfg.data['data'] + 'heatmaps/dataset.csv')
    with open(cfg.data['data'] + 'heatmaps/dataset.csv', "a", newline="") as f:
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=headers)
        if not file_exists:
            writer.writeheader()
        for row in rows:
            writer.writerow({'term': terms_names[row[0]], 'overlap': row[1], 'proteins': row[2], 'feature': feature})


def generate_heatmaps(feature):
    term_feature = pd.read_csv(cfg.data['data'] + 'heatmaps/dataset.csv')

    df_range = term_feature.loc[(term_feature.feature == feature)]

    # reshape term_feature dataset in proper format to create seaborn heatmap
    term_feature_df = df_range.pivot('term', 'overlap', 'proteins')
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=300)
    sns.set(font_scale=1)
    sns.heatmap(term_feature_df, cmap="BuPu", xticklabels=True, yticklabels=True)
    ax.set_ylabel('')
    # create seaborn heatmap
    plt.title(feature)
    plt.tick_params(labelsize=10)
    plt.savefig(cfg.data['data'] + 'heatmaps/range_' + feature + '.png', bbox_inches='tight')
    plt.close()


def draw_barplot(x, y):
    sns.barplot(x, y)
    plt.show()
