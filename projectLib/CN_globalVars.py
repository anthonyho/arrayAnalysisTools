# Anthony Ho, ahho@stanford.edu, 10/11/2015
# Last update 10/11/2016
"""Python module containing global variables for chemical nose project"""


import numpy as np
import pandas as pd
import seaborn as sns


# List of abbreviations of small molecules
sm_names = ['CA', 'GCA', 'TCA', 
            'DCA', 'GDCA', 'TDCA', 
            'CDCA', 'GCDCA', 
            'E1S', 'E2S', 'E3S', 
            'P5S', 'DHEAS', 
            'DHEAG', 'E2G'
            ]

# Number of small molecule
numSM = len(sm_names)

# Dictionary of directory number corresponding to each small molecule
sm_num = {'CA': 6, 
          'GCA': 7, 
          'TCA': 8, 
          'DCA': 9, 
          'GDCA': 14, 
          'TDCA': 15, 
          'CDCA': 10, 
          'GCDCA': 12, 
          'E1S': 2, 
          'E2S': 4, 
          'E3S': 3,
          'P5S': 1, 
          'DHEAS': 13, 
          'DHEAG': 11, 
          'E2G': 5
          }

# Dictionary of markers to be used for plotting for each small molecule
markers = {'CA': "^", 
           'GCA': "^", 
           'TCA': "^", 
           'DCA': "^", 
           'GDCA': "^", 
           'TDCA': "^", 
           'CDCA': "^", 
           'GCDCA': "^", 
           'E1S': "o", 
           'E2S': "o", 
           'E3S': "o",
           'P5S': "s", 
           'DHEAS': "s", 
           'DHEAG': "s", 
           'E2G': "o"
           }

# Dictionary of colors to be used for plotting for each small molecule
scolors = {'CA': sns.light_palette("purple", 5)[0], 
           'GCA': sns.light_palette("purple", 5)[1], 
           'TCA': sns.light_palette("purple", 5)[2], 
           'DCA': sns.light_palette("lightblue", 3)[0], 
           'GDCA': sns.light_palette("lightblue", 3)[1], 
           'TDCA': sns.light_palette("lightblue", 3)[2], 
           'CDCA': sns.dark_palette("blue", 5)[4], 
           'GCDCA': sns.dark_palette("blue", 5)[2], 
           'E1S': sns.light_palette("darkgreen", 6)[0], 
           'E2S': sns.light_palette("darkgreen", 6)[1], 
           'E3S': sns.light_palette("darkgreen", 6)[2],
           'P5S': sns.light_palette("darkgreen", 6)[3], 
           'DHEAS': sns.light_palette("darkgreen", 6)[4],
           'DHEAG': sns.light_palette("red", 4)[1], 
           'E2G': sns.light_palette("red", 4)[3]
           }

# Dictionary of category colors to be used for plotting for each small moleculew
catcolors = {'CA': sns.light_palette("purple", 5)[2], 
             'GCA': sns.light_palette("purple", 5)[2], 
             'TCA': sns.light_palette("purple", 5)[2], 
             'DCA': sns.light_palette("lightblue", 3)[2], 
             'GDCA': sns.light_palette("lightblue", 3)[2], 
             'TDCA': sns.light_palette("lightblue", 3)[2], 
             'CDCA': sns.dark_palette("blue", 5)[4], 
             'GCDCA': sns.dark_palette("blue", 5)[4], 
             'E1S': sns.light_palette("darkgreen", 6)[3], 
             'E2S': sns.light_palette("darkgreen", 6)[3], 
             'E3S': sns.light_palette("darkgreen", 6)[3],
             'P5S': sns.light_palette("darkgreen", 6)[3], 
             'DHEAS': sns.light_palette("darkgreen", 6)[3], 
             'DHEAG': sns.light_palette("red", 4)[3], 
             'E2G': sns.light_palette("red", 4)[3]
             }

# List of colors to be used for plotting for each small molecule
catcolors_list = [catcolors[currSM] for currSM in sm_names]

# Define class labels and categories
y1 = [0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 3]
catLabels1 = ['-cholic acids', '-deoxycholic acids', '-chenodeoxycholic acids', 'Estrogens', 'Androgens']
y2 = [0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4, 4]
catLabels2 = ['-cholic acids', '-deoxycholic acids', '-chenodeoxycholic acids', '3-sulfates', '3-glucuronides']
y3 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 1]
catLabels3 = ['Cholic acids', 'Estrogens', 'Androgens']
y4 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2]
catLabels4 = ['Cholic acids', '3-sulfates', '3-glucuronides']

# Define the highest concentration used (uM) in the titration series of sm
highest_conc = {'CA': 4096, 
                'GCA': 4096, 
                'TCA': 4096, 
                'DCA': 4096, 
                'GDCA': 4096, 
                'TDCA': 4096, 
                'CDCA': 2048, 
                'GCDCA': 4096, 
                'E1S': 1024, 
                'E2S': 1024, 
                'E3S': 1024,
                'P5S': 256, 
                'DHEAS': 1024, 
                'DHEAG': 512, 
                'E2G': 1024
                }

# Function to return a confusion matrix of true concentration of sm
def generate_true_conc_sm(conc, listSM=sm_names):
    return pd.DataFrame(np.identity(len(listSM)) * float(conc), index=listSM, columns=listSM)
    
# Define complex mixture id's
cm_id = ['cm_'+str(cm) for cm in range(1, 19)]

# Define true concentrations (uM) in complex mixtures
true_conc_cm = pd.DataFrame({'cm_1': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 44, 0, 235], 
                             'cm_2': [0, 0, 0, 0, 0, 0, 0, 0, 0, 228, 175, 105, 0, 0, 0], 
                             'cm_3': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 201, 301], 
                             'cm_4': [0, 0, 0, 0, 0, 0, 0, 0, 35, 0, 0, 0, 141, 0, 43], 
                             'cm_9':  [0, 0, 50, 0, 500, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             'cm_10': [15, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0, 0, 0, 0],
                             'cm_11': [0, 50, 0, 150, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0],
                             'cm_12': [1500, 0, 0, 0, 1200, 0, 0, 20, 0, 0, 0, 0, 0, 0, 0],
                             'cm_17': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 100, 0, 0], 
                             'cm_18': [0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 100]
                             }, 
                            index=sm_names)

