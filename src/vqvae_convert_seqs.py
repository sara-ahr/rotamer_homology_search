import pandas as pd
import numpy as np

import torch
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from math import sin, cos, isnan

import os



def predict(model, x):
    model.eval()
    with torch.no_grad():
        return model(torch.tensor(x, dtype=torch.float32)).detach().numpy()
    
def distance_matrix(a, b):
    return np.sqrt(np.sum((np.transpose(a) - b)**2, axis = 1))

def discretize(x, encoder, centroids):
    z = predict(encoder, x)
    return np.argmin(distance_matrix(z, centroids), axis = 0)

def insert_string(string, index, char):
    return string[:index] + char + string[index:]

AAs = ['ARG','ASN','ASP','CYS','GLU', 'GLN','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
concode = dict()
AAnum = 65
for a,i in enumerate(AAs):
    for x in range(3):
        if len(chr(AAnum)) == 1 and chr(AAnum) != '\\':
            concode[f'{i}{x}'] = chr(AAnum)
        else:    
            while len(chr(AAnum)) > 1 or chr(AAnum) == '\\':
                AAnum += 1
            concode[f'{i}{x}'] = chr(AAnum)
        AAnum += 1
if len(chr(AAnum)) == 1:
    concode['UNK'] = chr(AAnum)
else:    
    while len(chr(AAnum)) > 1:
        AAnum += 1
    concode['UNK'] = chr(AAnum)

for a in ['GLY','ALA']:
    AAnum += 1
    if len(chr(AAnum)) == 1:
        concode[a] = chr(AAnum)
    else:
        while len(chr(AAnum)) > 1:
            AAnum += AAnum
        concode[a] = chr(AAnum)



model = torch.load('/home/sarah/Downloads/Codes/encoder.pt')
states = np.loadtxt('/home/sarah/Downloads/Codes/states.txt')

integer_encoded = LabelEncoder().fit_transform(AAs)
integer_encoded = integer_encoded.reshape(len(integer_encoded),1)
onehot_encoded = OneHotEncoder(sparse_output=False).fit_transform(integer_encoded)
amino_np_code = dict(zip(AAs,onehot_encoded))


for file in os.listdir('/home/sarah/Downloads/pdb_example'):
    angles = pd.read_csv(os.path.join('/home/sarah/Downloads/pdb_example/', file))
    residf = angles.iloc[:,3]
    org_n = angles.shape[0]
    resils = []
    drop = []
    for i in range(residf.shape[0]):
        if residf.iloc[i] != 'UNK':
            resils.append(amino_np_code[residf.iloc[i]])
        else:
            drop.append(i)
    angles = angles.drop(drop)
    resinp = np.array(resils)

    tor = angles.iloc[:,5:]
    tor_ls = []
    for col in tor:
        sinls = []
        cosls = []
        for i in tor[col]:
            if isnan(i) != True:
                sinls.append(sin(i/180 * 3.1415))
                cosls.append(cos(i/180 * 3.1415))
            else:
                sinls.append(2)
                cosls.append(2)
        col = np.array([sinls, cosls])
        col = np.transpose(col)
        tor_ls.append(col)
    tor_np = np.hstack(tor_ls)
    all_inc_np = np.hstack([resinp, tor_np])
    output = np.apply_along_axis(discretize, 1, all_inc_np, model, states)
    output_char = ''.join([chr(i+65) for i in output])

    output_char = ''.join(concode)
    print(output_char)
    with open(f'/home/sarah/Downloads/pdb_con_seq_encoder_sep/{file[7:14]}.txt','w') as f:
        f.write(output_char)

    
    