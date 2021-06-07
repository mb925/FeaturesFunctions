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
    #     parse_file_protein_terms(term)

    for feature in features:
        for term in terms:
            print(term)
            # parse_file_proteins_features(term, feature)
            # generate_heatmaps_files(term, feature)
        generate_heatmaps(feature)


def parse_file_protein_terms(term):
    # print(term)

    with open(cfg.data['data'] + 'disprot_regions.json') as json_file:
        data = json.load(json_file)
        count = 0
        flag_term = False

        for el in data['data']:
            count += 1
            sequence_length = len(el['sequence'])
            term_column = [0] * sequence_length


            for region in el['regions']:
                if region['term_id'] == term:
                    flag_term = True
                    for i in range(region['start'] - 1, region['end'] -1):
                        term_column[i] = 1

            if flag_term:
                term_dir = cfg.data['data'] + 'proteins_terms_features/' + term + '/'
                if not os.path.exists(term_dir):
                    os.makedirs(term_dir)
                df = pd.DataFrame(term_column, columns=['term'])
                df.to_csv(term_dir + el['acc'] + '.tsv', sep='\t')

        if flag_term == False:
            print(term + ' not found')

def parse_file_proteins_features(term, feature):
    count = 0
    df_features = pd.read_csv(cfg.data['data'] + 'mobidb_features.tsv', sep='\t')
    print(feature)

    term_dir = cfg.data['data'] + 'proteins_terms_features/' + term + '/'
    for filename in os.listdir(term_dir):
        count += 1
        acc = filename.split('.')[0]
        regions = df_features.loc[(df_features.acc == acc) & (df_features.feature == feature)]['start..end'].values
        if len(regions) > 0:
            regions = regions[0].split(',')

        df_file = pd.read_csv(term_dir + filename, sep='\t', index_col=[0])

        # will be zero if feature is missing or if acc is missing in the mobidb file
        feature_column = [0] * len(df_file)
        for region in regions:
            start = region.split('..')[0]
            end = region.split('..')[1]
            for i in range(int(start) - 1, int(end) - 1):
                feature_column[i] = 1
        df_file[feature] = feature_column

        df_file.to_csv(term_dir + filename, sep='\t')

def calculate_jaccard_index(term, feature):

    dict_jaccard = {term + '-' + feature: []}
    term_dir = cfg.data['data'] + 'proteins_terms_features/' + term + '/'

    for filename in os.listdir(term_dir):
        print(filename)
        df_file = pd.read_csv(term_dir + filename, sep='\t', index_col=[0])
        intersection = df_file.loc[(df_file.term == 1) & (df_file[feature] == 1)]
        union = df_file.loc[(df_file.term == 1) | (df_file[feature] == 1)]
        jaccard = 0
        if len(intersection) != 0:
            jaccard = round(len(intersection) / len(union), 2)
            # print(jaccard)
        if jaccard != 0:
            dict_jaccard[term + '-' + feature].append(jaccard)


    y = np.array(dict_jaccard[term + '-' + feature])
    plt.hist(y)
    ax = plt.gca()  # get axis handle
    p = ax.patches
    heights = [patch.get_height() for patch in p]
    norm = np.linalg.norm(heights)

    if norm == 0:
        normal_heights = [0]*10
    else:
        normal_heights = heights / norm


    plt.close()
    row = []
    for j in range(10):
        row.append(term)
        row.append(j + 1)
        row.append(round(normal_heights[j], 2))

    row_splitted = [row[i:i + 3] for i in range(0, len(row), 3)]
    return row_splitted

def generate_heatmaps_files(term, feature):
    rows = calculate_jaccard_index(term, feature)
    headers = ['term', 'jaccard', 'proteins', 'feature']
    file_exists = os.path.isfile(cfg.data['data'] + 'heatmaps/dataset.csv')
    with open(cfg.data['data'] + 'heatmaps/dataset.csv', "a", newline="") as f:
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=headers)
        if not file_exists:
            writer.writeheader()
        for row in rows:
            writer.writerow({'term': terms_names[row[0]], 'jaccard': row[1], 'proteins': row[2], 'feature': feature})

def generate_heatmaps(feature):
    term_feature = pd.read_csv(cfg.data['data'] + 'heatmaps/dataset.csv')


    df_range = term_feature.loc[(term_feature.feature == feature)]

    # reshape term_feature dataset in proper format to create seaborn heatmap
    term_feature_df = df_range.pivot('term', 'jaccard', 'proteins')
    print(term_feature_df)

    sns.heatmap(term_feature_df, cmap="BuPu", xticklabels=True, yticklabels=True)

  # create seaborn heatmap
    plt.title(feature)
    plt.tick_params(labelsize=5)
    plt.savefig(cfg.data['data'] + 'heatmaps/range_' + feature + '.png', bbox_inches='tight')
    plt.close()



def draw_barplot(x, y):
    sns.barplot(x, y)
    plt.show()
