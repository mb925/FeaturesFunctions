import numpy as np
import matplotlib.pyplot as plt
import parser
import obo_parser
import config as cfg



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    obo_file = cfg.data['data'] + 'idpontology_disprot.obo'

    ontology_tree = obo_parser.parse_obo(obo_file)
    # list_terms = parser.visit_tree(ontology_tree, 3)
    # print(list_terms)
    # list_features = ['prediction-negative_polyelectrolyte-mobidb_lite_sub', 'prediction-polar-mobidb_lite_sub']
    # # list_features = ['prediction-low_complexity-merge']
    # parser.create_protein_terms_features_table(list_terms, list_features, 'level_3')
    # print(list_features)

    # level = 1
    # list_terms = parser.visit_tree(ontology_tree, 3) + parser.visit_tree(ontology_tree, 2) + parser.visit_tree(ontology_tree, 1)

    level = 2
    relationships = obo_parser.get_relationships(parser.visit_tree(ontology_tree, 2), ontology_tree)
    list_terms = parser.visit_tree(ontology_tree, 3) + parser.visit_tree(ontology_tree, 2)
    # print(list_terms)
    # list_terms = parser.visit_tree(ontology_tree, 3) + parser.visit_tree(ontology_tree, 2)
    # list_terms = parser.visit_tree(ontology_tree, 2)

    # level = 3
    # list_terms = parser.visit_tree(ontology_tree, 3)



    list_features = ['curated-disorder-disprot' ]
    parser.create_protein_terms_features_table(list_terms, list_features, 'level_' + str(level), relationships)
    # print(list_features)
