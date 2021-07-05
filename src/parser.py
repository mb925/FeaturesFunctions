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

    '00078': {'name':'pre-molten globule', 'quantity':''},
    '00077': {'name':'molten globule', 'quantity':''},
    '00051': {'name':'disorder to molten globule', 'quantity':''},
    '00052': {'name':'disorder to pre-molten globule', 'quantity':''},
    '00053': {'name':'molten globule to order', 'quantity':''},
    '00055': {'name':'pre-molten globule to order', 'quantity':''},
    '00058': {'name':'Order to pre-molten globule', 'quantity':''},
    '00060': {'name':'Pre-molten globule to molten globule', 'quantity':''},
    '00063': {'name':'Protein binding', 'quantity':''},
    '00065': {'name':'DNA binding', 'quantity':''},
    '00066': {'name':'RNA binding', 'quantity':''},
    '00071': {'name':'Lipid binding', 'quantity':''},
    '00073': {'name':'Metal ion binding', 'quantity':''},
    '00074': {'name':'Small molecule binding', 'quantity':''},
    '00002': {'name':'Flexible linker/spacer', 'quantity':''},
    '00003': {'name':'Entropic bristle', 'quantity':''},
    '00005': {'name':'Entropic spring', 'quantity':''},
    '00006': {'name':'Structural mortar', 'quantity':''},
    '00007': {'name':'Self-transport through channel', 'quantity':''},
    '00009': {'name':'Assembler', 'quantity':''},
    '00011': {'name':'Targeting', 'quantity':''},
    '00012': {'name':'Tethering', 'quantity':''},
    '00014': {'name':'Neutralization of toxic molecules', 'quantity':''},
    '00015': {'name':'Metal binding', 'quantity':''},
    '00016': {'name':'Water storage', 'quantity':''},
    '00019': {'name':'Cis-regulatory element', 'quantity':''},
    '00020': {'name':'Disassembler', 'quantity':''},
    '00021': {'name':'Activator', 'quantity':''},
    '00022': {'name':'DNA bending', 'quantity':''},
    '00023': {'name':'DNA unwinding', 'quantity':''},
    '00025': {'name':'Regulation of phosphorylation', 'quantity':''},
    '00026': {'name':'Regulation of acetylation', 'quantity':''},
    '00027': {'name':'Regulation of methylation', 'quantity':''},
    '00028': {'name':'Regulation of glycosylation', 'quantity':''},
    '00029': {'name':'Regulation of ubiquitination', 'quantity':''},
    '00032': {'name':'Regulation of palmitoylation', 'quantity':''},
    '00033': {'name':'Regulation of limited proteolysis', 'quantity':''},
    '00034': {'name':'Regulation of ADP-ribosylation', 'quantity':''},
    '00037': {'name':'Space filling', 'quantity':''},
    '00038': {'name':'Entropic exclusion', 'quantity':''},
    '00041': {'name':'Liquid-liquid phase separation', 'quantity':''},
    '00042': {'name':'Hydrogel', 'quantity':''},
    '00043': {'name':'Aggregate', 'quantity':''},
    '00046': {'name':'Amyloid', 'quantity':''},
    '00047': {'name':'Prion', 'quantity':''}

}
dict_tree = {}
dict_tree_set = {}

def access(tree, value):

    while len(tree[value]['parents']) > 0:
        print('my dad!')
        dict_tree[value].append(tree[value]['parents'][0])

        value = tree[value]['parents'][0]
        print(value)
        access(tree, value)


def obtain_terms(tree):

    for value in tree:
        dict_tree[value] = []
        print('I am: ' + value)
        access(tree, value)

    print(dict_tree)



def build_dictionary_ontology():
    file = cfg.data['data'] + 'proteins_terms_features/terms_features.tsv'
    df_terms_features = pd.read_csv(file, sep='\t')
    for term in terms_names:
        file_count = df_terms_features.loc[(df_terms_features.term_code == terms_names[term]['name'])].acc.unique()
        terms_names[term]['quantity'] = str(len(file_count))
    print(terms_names)


def create_protein_terms_features_table(terms, features):
    # for term in terms:
    #     parse_file_protein_terms(term)
    #
    # build_dictionary_ontology()

    for feature in features:
        # add features to files
        # parse_file_proteins_features(terms, feature)
        # generate_heatmap_dataset(terms, feature)
        generate_heatmaps(feature)
    # # I could erase all generated files with a command here at the end


def parse_file_protein_terms(term):
    print(term)

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
                    'term_code': terms_names[term]['name'],
                    'term': dict_terms[el['acc']]
                })
                df = df.append(df_el)

        if the_term_is_in_the_whole_file == False:
            print(term + ' not found')

        file = cfg.data['data'] + 'proteins_terms_features/terms_features.tsv'
        df.to_csv(file, mode='a',  header=(not os.path.exists(file)), sep='\t', index=False)




def parse_file_proteins_features(terms, feature):
    full_feature_column = []
    terms_features_file = cfg.data['data'] + 'proteins_terms_features/terms_features.tsv'

    df_terms_features = pd.read_csv(terms_features_file, sep='\t')

    for term in terms:

        df_features = pd.read_csv(cfg.data['data'] + 'mobidb_features.tsv', sep='\t')
        print(feature)

        # calculate features for each acc
        accessions = df_terms_features.loc[(df_terms_features.term_code == terms_names[term]['name'])].acc.unique()
        for acc in accessions:

            sequence_length = df_terms_features.loc[(df_terms_features.acc == acc) & (df_terms_features.term_code == terms_names[term]['name']) ]
            sequence_length = len(sequence_length.index)

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
    terms_features_file = cfg.data['data'] + 'proteins_terms_features/terms_features.tsv'
    df_terms_features = pd.read_csv(cfg.data['data'] + 'proteins_terms_features/terms_features.tsv', sep='\t')

    accessions = df_terms_features.loc[(df_terms_features.term_code == terms_names[term]['name'])].acc.unique()
    for acc in accessions:
        df_file = pd.read_csv(terms_features_file, sep='\t', index_col=[0])
        intersection = df_file.loc[(df_file.acc == acc) & (df_file.term == 1) & (df_file[feature] == 1)]
        union = df_file.loc[(df_file.acc == acc) & ((df_file.term == 1) | (df_file[feature] == 1))]
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
        row.append(term)
        row.append(x + 1)
        # normalizing the occurrences to the range 0 1
        row.append(total_occurrences[x] / sum(total_occurrences))

    row_splitted = [row[i:i + 3] for i in range(0, len(row), 3)]

    return row_splitted


def generate_heatmap_dataset(terms, feature):
    for term in terms:
        rows = calculate_jaccard_index(term, feature)
        headers = ['term', 'overlap', 'proteins', 'feature']
        file_exists = os.path.isfile(cfg.data['data'] + 'heatmaps/dataset.csv')
        with open(cfg.data['data'] + 'heatmaps/dataset.csv', "a", newline="") as f:
            writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=headers)
            if not file_exists:
                writer.writeheader()
            for row in rows:
                writer.writerow({'term': terms_names[row[0]]['name'] + ' (' + terms_names[row[0]]['quantity'] + ')', 'overlap': row[1], 'proteins': row[2], 'feature': feature}) #  todo: add quantity of annotation


def generate_heatmaps(feature):
    term_feature = pd.read_csv(cfg.data['data'] + 'heatmaps/dataset.csv')

    df_range = term_feature.loc[(term_feature.feature == feature)]

    # reshape term_feature dataset in proper format to create seaborn heatmap
    term_feature_df = df_range.pivot('term', 'overlap', 'proteins')
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=300)
    sns.set(font_scale=1)
    sns.heatmap(term_feature_df, cmap="BuPu", xticklabels=True, yticklabels=True, vmin=0, vmax=1)
    ax.set_ylabel('')
    # create seaborn heatmap
    plt.title(feature)
    plt.tick_params(labelsize=10)
    plt.savefig(cfg.data['data'] + 'heatmaps/range_' + feature + '.png', bbox_inches='tight')
    plt.close()


def draw_barplot(x, y):
    sns.barplot(x, y)
    plt.show()
