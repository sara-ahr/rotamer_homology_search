import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
import os

from math import sin, cos, isnan

os.chdir('/home/sarah/Downloads/Codes')
aligned_angles = pd.read_csv('alignment_cut_list.csv')

amino_acid_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

integer_encoded = LabelEncoder().fit_transform(amino_acid_codes)
integer_encoded = integer_encoded.reshape(len(integer_encoded),1)
onehot_encoded = OneHotEncoder(sparse_output=False).fit_transform(integer_encoded)
amino_np_code = dict(zip(amino_acid_codes,onehot_encoded))


print(aligned_angles.shape)
resi1resi2 = aligned_angles.iloc[:,[1,9]]
resi1ls = []
resi2ls = []
drop = []
for i in range(resi1resi2.shape[0]):
    print(i)
    if resi1resi2.iloc[i,0] != 'UNK' and resi1resi2.iloc[i,1] != 'UNK':
        resi1ls.append(amino_np_code[resi1resi2.iloc[i,0]])
        resi2ls.append(amino_np_code[resi1resi2.iloc[i,1]])
    else:
        drop.append(i)
print(drop)
print(len(resi1ls))
aligned_angles = aligned_angles.drop(drop)
print(aligned_angles.shape)

resi1np = np.array(resi1ls)
print(resi1np.shape)
resi2np = np.array(resi2ls)
print(resi2np.shape)


tor_1ls = []
tor_1 = aligned_angles.iloc[:,2]
print(tor_1)
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
tor_2 = aligned_angles.iloc[:,10]
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

all_inc_np = np.hstack([resi1np, tor_1np, resi2np, tor_2np])
np.save('alignment_modified_vectors_c1.npy',all_inc_np)




