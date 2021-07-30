import parser
import obo_parser
import config as cfg


if __name__ == '__main__':
    obo_file = cfg.data['data'] + 'idpontology_disprot.obo'
    # list_features = ['prediction-low_complexity-merge', 'prediction-lip-anchor', 'prediction-proline_rich-mobidb_lite_sub',
    #                  'prediction-negative_polyelectrolyte-mobidb_lite_sub', 'prediction-positive_polyelectrolyte-mobidb_lite_sub',
    #                  'prediction-polyampholyte-mobidb_lite_sub', 'prediction-polar-mobidb_lite_sub']
    # list_features = ['prediction-low_complexity-merge']
    list_features = ['curated-disorder-disprot']
    # generate the tree
    ontology_tree = obo_parser.parse_obo(obo_file)

    # execute the level/s you need

    level = 1
    relationships = obo_parser.get_relationships(parser.visit_tree(ontology_tree, 1) + parser.visit_tree(ontology_tree, 2), ontology_tree)
    list_terms = parser.visit_tree(ontology_tree, 3) + parser.visit_tree(ontology_tree, 2) + parser.visit_tree(ontology_tree, 1)
    parser.create_protein_terms_features_table(list_terms, list_features, 'level_' + str(level), relationships)

    level = 2
    relationships = obo_parser.get_relationships(parser.visit_tree(ontology_tree, 2), ontology_tree)
    list_terms = parser.visit_tree(ontology_tree, 3) + parser.visit_tree(ontology_tree, 2)
    parser.create_protein_terms_features_table(list_terms, list_features, 'level_' + str(level), relationships)


    level = 3
    list_terms = parser.visit_tree(ontology_tree, 3)
    parser.create_protein_terms_features_table(list_terms, list_features, 'level_' + str(level), [])

    # print(list_features)

    # if you run code to generate datasets for all 3 levels, you can create heatmaps
    for feature in list_features:
        parser.generate_heatmaps(feature)
