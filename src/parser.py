import csv
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
    '00001': {'name':'Entropic chain', 'quantity':''},
    '00002': {'name':'Flexible linker/spacer', 'quantity':''},
    '00003': {'name':'Entropic bristle', 'quantity':''},
    '00004': {'name':'Entropic clock', 'quantity':''},
    '00005': {'name':'Entropic spring', 'quantity':''},
    '00006': {'name':'Structural mortar', 'quantity':''},
    '00007': {'name':'Self-transport through channel', 'quantity':''},
    '00008': {'name':'Molecular recognition assembler', 'quantity':''},
    '00009': {'name':'Assembler', 'quantity':''},
    '00010': {'name':'Localization', 'quantity':''},
    '00011': {'name':'Targeting', 'quantity':''},
    '00012': {'name':'Tethering', 'quantity':''},
    '00013': {'name':'Molecular recognition scavenger', 'quantity':''},
    '00014': {'name':'Neutralization of toxic molecules', 'quantity':''},
    '00015': {'name':'Metal binding', 'quantity':''},
    '00016': {'name':'Water storage', 'quantity':''},
    '00017': {'name':'Molecular recognition effector', 'quantity':''},
    '00018': {'name':'Inhibitor', 'quantity':''},
    '00019': {'name':'Cis-regulatory element', 'quantity':''},
    '00020': {'name':'Disassembler', 'quantity':''},
    '00021': {'name':'Activator', 'quantity':''},
    '00022': {'name':'DNA bending', 'quantity':''},
    '00023': {'name':'DNA unwinding', 'quantity':''},
    '00024': {'name':'Molecular recognition display site', 'quantity':''},
    '00025': {'name':'Regulation of phosphorylation', 'quantity':''},
    '00026': {'name':'Regulation of acetylation', 'quantity':''},
    '00027': {'name':'Regulation of methylation', 'quantity':''},
    '00028': {'name':'Regulation of glycosylation', 'quantity':''},
    '00029': {'name':'Regulation of ubiquitination', 'quantity':''},
    '00030': {'name':'Regulation of fatty acylation', 'quantity':''},
    '00031': {'name':'Regulation of myristoylation', 'quantity':''},
    '00032': {'name':'Regulation of palmitoylation', 'quantity':''},
    '00033': {'name':'Regulation of limited proteolysis', 'quantity':''},
    '00034': {'name':'Regulation of ADP-ribosylation', 'quantity':''},
    '00035': {'name':'Molecular recognition chaperone', 'quantity':''},
    '00036': {'name':'Protein detergent', 'quantity':''},
    '00037': {'name':'Space filling', 'quantity':''},
    '00038': {'name':'Entropic exclusion', 'quantity':''},
    '00039': {'name':'Entropy transfer', 'quantity':''},
    '00040': {'name':'Biological condensation', 'quantity':''},
    '00041': {'name':'Liquid-liquid phase separation', 'quantity':''},
    '00042': {'name':'Hydrogel', 'quantity':''},
    '00043': {'name':'Aggregate', 'quantity':''},
    '00044': {'name':'Granule', 'quantity':''},
    '00045': {'name':'Cellular puncta', 'quantity':''},
    '00046': {'name':'Amyloid', 'quantity':''},
    '00047': {'name':'Prion', 'quantity':''},
    '00048': {'name':'Signaling complex', 'quantity':''}

}

def create_protein_terms_features_table(terms, features, level, relationships):
    read_disprot_terms()
    # for term in terms:
    #     terms = parse_file_protein_terms(terms, term, level)

    build_dictionary_ontology(level) #needed to set quantity

    for feature in features:
    #     # add features to files
    #     parse_file_proteins_features(terms, feature, level)
        generate_heatmap_dataset(terms, feature, level, relationships)
        # generate_heatmaps(feature, level)
    # I could erase all generated files with a command here at the end

def visit_tree(tree, target):
    tree_visit = []
    levels_tree = {}
    registered = [] # values that are registered into the levels_tree

    key, value = 'DO:00000', tree['DO:00000'] # choose which branch to keep
    level = 0
    node = [key, value, level]
    tree_visit.append(node)
    levels_tree[node[2]] = [key.split(':')[0]]
    registered.append(key)
    while len(tree_visit) > 0:
        node = tree_visit.pop()
        if len(node[1]['parents']) > 0:
            parent_key = node[1]['parents'][0]
            if parent_key not in registered:
                    parent_value = tree[parent_key]
                    parent_level = node[2] - 1
                    if parent_level in levels_tree:
                        levels_tree[parent_level].append(parent_key.split(':')[1])
                    else:
                        levels_tree[parent_level] = [parent_key.split(':')[1]]

                    registered.append(parent_key)
                    tree_visit.append([parent_key, parent_value, parent_level])

        for child in node[1]['children']:
            if child not in registered:

                child_value = tree[child]
                child_level = node[2] + 1
                if child_level in levels_tree:
                    levels_tree[child_level].append(child.split(':')[1])
                else:
                    levels_tree[child_level] = [child.split(':')[1]]

                registered.append(child)
                tree_visit.append([child, child_value, child_level])

    list_terms = levels_tree[target]
    return list_terms



def build_dictionary_ontology(level):
    file = cfg.data['data'] + 'proteins_terms_features/' + level + '/terms_features.tsv'
    df_terms_features = pd.read_csv(file, sep='\t')
    for term in terms_names:
        file_count = df_terms_features.loc[(df_terms_features.term_code == terms_names[term]['name'])].acc.unique()
        terms_names[term]['quantity'] = str(len(file_count)) # todo add child when neeeded

def read_disprot_terms():
    file = cfg.data['data'] + '/disprot_terms.tsv'
    df = pd.read_csv(file, sep='\t')
    df = df.loc[(df.acc == 'P03265')]
    # print(df)

    # Plot
    # fig, axes = plt.subplots(1, 5, figsize=(10,2.5), dpi=100, sharex=True, sharey=True)
    # colors = ['tab:red', 'tab:blue', 'tab:green', 'tab:pink', 'tab:olive']
    #
    # for i, (ax, cut) in enumerate(zip(axes.flatten(), df.cut.unique())):
    #     x = df.loc[df.cut==cut, 'depth']
    #     ax.hist(x, alpha=0.5, bins=100, density=True, stacked=True, label=str(cut), color=colors[i])
    #     ax.set_title(cut)
    #
    # plt.suptitle('Probability Histogram of Diamond Depths', y=1.05, size=16)
    # ax.set_xlim(50, 70); ax.set_ylim(0, 1);
    # plt.tight_layout()




def parse_file_protein_terms(list_terms, term, level):
    # print(term)

    with open(cfg.data['data'] + 'disprot_regions.json') as json_file:
        data = json.load(json_file)
        dict_terms = {}
        df = pd.DataFrame()
        the_term_is_in_the_whole_file = False

        for el in data['data']:
            the_term_is_in_the_current_row = False

            sequence_length = len(el['sequence'])

            for region in el['regions']:
                if region['term_namespace'] == 'Disorder function':
                    if region['term_id'] == term:
                        the_term_is_in_the_whole_file = True
                        the_term_is_in_the_current_row = True

                        for i in range(region['start'] - 2, region['end']):
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
            list_terms.remove(term)
        else:
            file = cfg.data['data'] + 'proteins_terms_features/' + level + '/terms_features.tsv' # todo: create dir
            df.to_csv(file, mode='a',  header=(not os.path.exists(file)), sep='\t', index=False)
        return list_terms

def parse_file_proteins_features(terms, feature, level):
    full_feature_column = []
    terms_features_file = cfg.data['data'] + 'proteins_terms_features/' + level + '/terms_features.tsv' # todo: create dir

    df_terms_features = pd.read_csv(terms_features_file, sep='\t')

    for term in terms:

        df_features = pd.read_csv(cfg.data['data'] + 'mobidb_features.tsv', sep='\t')
        # print(feature)

        # calculate features for each acc
        # print(terms_names)
        accessions = df_terms_features.loc[(df_terms_features.term_code == terms_names[term]['name'])].acc.unique()
        for acc in accessions:

            sequence_length = df_terms_features.loc[(df_terms_features.acc == acc) & (df_terms_features.term_code == terms_names[term]['name']) ]
            sequence_length = len(sequence_length.index)

            regions = df_features.loc[(df_features.acc == acc) & (df_features.feature == feature)]['start..end'].values
            if len(regions) > 0:
                regions = regions[0].split(',')
                # print(regions)

            # write features to file
            # will be zero if feature is missing or if acc is missing in the mobidb file
            feature_column = [0] * sequence_length
            for region in regions:
                start = region.split('..')[0]
                end = region.split('..')[1]
                for i in range(int(start) - 2, int(end)):
                    # print(i)
                    feature_column[i] = 1
            full_feature_column = full_feature_column + feature_column

    df_terms_features[feature] = full_feature_column
    df_terms_features.to_csv(terms_features_file, sep='\t')


def calculate_jaccard_index(term, feature, level, relationship):

    dict_jaccard = {term + '-' + feature: []}
    terms_features_file = cfg.data['data'] + 'proteins_terms_features/' + level + '/terms_features.tsv'
    df_terms_features = pd.read_csv(cfg.data['data'] + 'proteins_terms_features/' + level + '/terms_features.tsv', sep='\t')
    jaccard_values = []
    # consider children as well
    # HERE IS WHERE FATHER NODES ENGLOBE CHILDREN
    if term in relationship:
        accessions = df_terms_features.loc[(df_terms_features.term_code == terms_names[term]['name'])].acc.unique().tolist()
        for child in relationship[term]:
            acc_child = df_terms_features.loc[(df_terms_features.term_code == terms_names[child]['name'])].acc.unique().tolist()
            accessions = accessions + acc_child
        print(accessions)
    else:
        accessions = df_terms_features.loc[(df_terms_features.term_code == terms_names[term]['name'])].acc.unique().tolist()
    left = 0
    for acc in accessions:
        # print('acc: ' + acc)
        df_file = pd.read_csv(terms_features_file, sep='\t', index_col=[0])
        intersection = df_file.loc[(df_file.acc == acc) & (df_file.term == 1) & (df_file[feature] == 1)]

        # union = df_file.loc[(df_file.acc == acc) & ( df_file[feature] == 1)] #feature
        union = df_file.loc[(df_file.acc == acc) & (df_file.term == 1)] #function
        # union = df_file.loc[(df_file.acc == acc) & ((df_file.term == 1) | (df_file[feature] == 1))]
        jaccard = 0
        if len(intersection) != 0:
            jaccard = round(len(intersection) / len(union), 2)
            # print(jaccard)
        if jaccard != 0:
            left += 1
            # print(jaccard)
        dict_jaccard[term + '-' + feature].append(round(jaccard, 1))
        jaccard_values.append(round(jaccard, 1))

    row = []
    # how many times each step happens
    total_occurrences = []

    for j in np.arange(0.1, 1.1, 0.1): # [0.1, 0.2, 0.3, 0.4 ...] -- skipping all 0 values
        jaccard_occurrence = dict_jaccard[term + '-' + feature].count(round(j, 1))
        total_occurrences.append(jaccard_occurrence) # number of proteins (accession) in that range
    terms_names[term]['left'] = str(left)

    for x in range(10):
        row.append(term)
        row.append(x + 1)
        # normalizing the occurrences to the range 0 1
        if sum(total_occurrences) == 0:
            row.append(0)
        else:
            row.append(total_occurrences[x] / sum(total_occurrences))

    row_splitted = [row[i:i + 3] for i in range(0, len(row), 3)]

    # print('length')
    # print(len(jaccard_values))
    return row_splitted, jaccard_values


def generate_heatmap_dataset(terms, feature, level, relationship):

    for term in terms:
        # print(term)
        # if you are a child I'll skip you
        # I'm going to consider u anyway through relationship with father when I calculate jaccard

        if term in sum(relationship.values(), []):
            continue
        rows, jaccard_values = calculate_jaccard_index(term, feature, level, relationship)
        headers = ['term', 'overlap', 'proteins', 'feature']
        file_exists = os.path.isfile(cfg.data['data'] + 'heatmaps/' + level + '/dataset.csv')
        with open(cfg.data['data'] + 'heatmaps/' + level + '/dataset.csv', "a", newline="") as f:
            writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=headers)
            if not file_exists:
                writer.writeheader()
            for row in rows:
                # print(terms_names)
                # print(row[0])
                if int(terms_names[row[0]]['quantity']) > 10: # todo: && child
                    writer.writerow({'term': terms_names[row[0]]['name'] + ' (' + terms_names[row[0]]['left'] + '/' + terms_names[row[0]]['quantity'] + ')', 'overlap': row[1], 'proteins': row[2], 'feature': feature})

        # if len(jaccard_values) > 0:
        #     dict_jac[term + '_' + feature] = jaccard_values

# def draw_histograms(total, fig, count, jaccard_values):
#     print(total)
#     print(count)
#     ax1 = fig.add_subplot(1, total, count)
#     ax1.set_ylabel('#proteins', fontsize=20)
#     ax1.set_xlabel('overlap', fontsize=20)
#     ax1.hist(jaccard_values, bins=20)
#     plt.show()

def generate_heatmaps(feature, level):
    term_feature = pd.read_csv(cfg.data['data'] + 'heatmaps/' + level + '/dataset.csv')

    df_range = term_feature.loc[(term_feature.feature == feature)]

    # reshape term_feature dataset in proper format to create seaborn heatmap
    term_feature_df = df_range.pivot('term', 'overlap', 'proteins')
    fig, ax = plt.subplots(1, 1, figsize=(8, 3), dpi=300)
    sns.set(font_scale=1)
    sns.heatmap(term_feature_df, cmap="BuPu", xticklabels=True, yticklabels=True, vmin=0, vmax=1)
    ax.set_ylabel('')
    # create seaborn heatmap
    plt.title(feature)
    plt.tick_params(labelsize=10)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()
    plt.savefig(cfg.data['data'] + 'heatmaps/' + level + '/range_' + feature + '.png', bbox_inches='tight')
    plt.close()


def draw_barplot(x, y):
    sns.barplot(x, y)
    plt.show()
