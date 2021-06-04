import os

import numpy as np
import pandas as pd
import config as cfg
import json
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

def create_protein_terms_features_table(terms, features):

    # for term in terms:
    #     parse_file_protein_terms(term)

    for feature in features:
        print(feature)
        for term in terms:
            print(term)
            # parse_file_proteins_features(term, feature)
            generate_heatmaps(term)


def parse_file_protein_terms(term):
    print(term)

    with open(cfg.data['data'] + 'disprot_regions.json') as json_file:
        data = json.load(json_file)
        count = 0
        for el in data['data']:
            count += 1
            sequence_length = len(el['sequence'])
            term_column = [0] * sequence_length
            flag_term = False


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
        df_file = pd.read_csv(term_dir + filename, sep='\t', index_col=[0])
        intersection = df_file.loc[(df_file.term == 1) & (df_file[feature] == 1)]
        union = df_file.loc[(df_file.term == 1) | (df_file[feature] == 1)]
        jaccard = round(len(intersection) / len(union), 2)
        # print(jaccard)
        if jaccard != 0:
            dict_jaccard[term + '-' + feature].append(jaccard)

    print(dict_jaccard[term + '-' + feature])
    y = np.array(dict_jaccard[term + '-' + feature])
    plt.hist(y)
    ax = plt.gca()  # get axis handle
    p = ax.patches
    heights = [patch.get_height() for patch in p]
    print(heights)
    plt.close()
    row = []
    for i in range(10):
        row.append(term)
        row.append(i + 1)
        row.append(heights[i])
    return row

def generate_heatmaps(term, feature):
    row = calculate_jaccard_index(term, feature)
    if os.path.isfile(feature + '.csv'):
        # scrivi il nuovo termine
        with open(cfg.data['data'] + 'heatmaps/' + feature + '.csv', 'a') as file:
            for line in row:
                file.write(line)
                file.write('\n')
    else:
        # crei file e scrivi il nuovo termine
        with open(cfg.data['data'] + 'heatmaps' + feature + '.csv', 'w') as file:
            for line in row:
                file.write(line)
                file.write('\n')



def draw_barplot(x, y):
    sns.barplot(x, y)
    plt.show()
