import numpy as np
import matplotlib.pyplot as plt
import parser
import seaborn as sns # for data visualization



# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    list_terms = ['00076']
    list_features = ['prediction-low_complexity-merge']
    # list_terms = ['00076', '00079']
    # list_features = ['prediction-low_complexity-merge', 'prediction-proline_rich-mobidb_lite_sub']
    parser.create_protein_terms_features_table(list_terms, list_features)

    # flight = sns.load_dataset('flights')  # load flights datset from GitHub seaborn repository
    #
    # # reshape flights dataset in proper format to create seaborn heatmap
    # flights_df = flight.pivot('month', 'year', 'passengers')
    #
    # sns.heatmap(flights_df)# create seaborn heatmap
    #
    # plt.show()
    # plt.close()
