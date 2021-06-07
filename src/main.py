import numpy as np
import matplotlib.pyplot as plt
import parser
import seaborn as sns # for data visualization



# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    list_terms = ['00077', '00078',  '00051',  '00052',  '00053', '00055',
                  '00058',  '00060', '00063', '00065',
                  '00066', '00071', '00073', '00074', '00002', '00003',
                  '00005', '00006', '00007', '00009', '00011', '00012',
                  '00014', '00015', '00016', '00019', '00020', '00021',
                  '00022', '00023', '00025', '00026', '00027', '00028',
                  '00029', '00032', '00033', '00034', '00037', '00038',
                  '00041', '00042', '00043', '00046', '00047']
    list_features = ['prediction-low_complexity-merge', 'prediction-proline_rich-mobidb_lite_sub',
                     'prediction-lip-anchor']
    parser.create_protein_terms_features_table(list_terms, list_features)


