import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from math import sin, cos, isnan

os.chdir('/home/sarah/Downloads/Codes')
aligned_angles = pd.read_csv('alignment_modified_list.csv')
os.chdir('/home/sarah/Downloads/Alignment_vectors_chi1chi2')



amino_acid_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

dict_list_vectors = dict(zip(amino_acid_codes, [[] for i in range(20)]))

for row in range(aligned_angles.shape[0]):
    if aligned_angles.iloc[row,1] == aligned_angles.iloc[row,9] and aligned_angles.iloc[row,1] != 'UNK':
        dict_list_vectors[aligned_angles.iloc[row,1]].append([i for i in aligned_angles.iloc[row,:]])

dict_dataframes = {}
for key in dict_list_vectors.keys():
    dict_dataframes[key] = pd.DataFrame(dict_list_vectors[key])

for key in dict_dataframes.keys():
    value = dict_dataframes[key]
    tor_1ls = []
    tor_1 = value.iloc[:,[2,4]]
    for col in tor_1:
        sinls = []
        cosls = []
        for i in tor_1[col]:
            if isnan(i) != True:
                sinls.append(sin(i/180 * 3.1415))
                cosls.append(cos(i/180 * 3.1415))
            else:
                sinls.append(2)
                cosls.append(2)
        col = np.array([sinls, cosls])
        col = np.transpose(col)
        tor_1ls.append(col)
    tor_1np = np.hstack(tor_1ls)      


    tor_2ls = []
    tor_2 = value.iloc[:,[10,12]]
    for col in tor_2:
        sinls = []
        cosls = []
        for i in tor_2[col]:
            if isnan(i) != True:
                sinls.append(sin(i/180 * 3.1415))
                cosls.append(cos(i/180 * 3.1415))
            else:
                sinls.append(2)
                cosls.append(2)
        col = np.array([sinls, cosls])
        col = np.transpose(col)
        tor_2ls.append(col)
    tor_2np = np.hstack(tor_2ls)

    all_inc_np = np.hstack([tor_1np, tor_2np])
    print(all_inc_np)
    np.save(f'alignment_modified_vectors_{key}.npy',all_inc_np)
