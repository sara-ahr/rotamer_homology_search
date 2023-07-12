import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from math import sin, cos, isnan
import torch


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


models = {}
for model in os.listdir('/home/sarah/Downloads/AA_encoders/Encoders'):
    models[model[-6:-3]] = torch.load(os.path.join('/home/sarah/Downloads/AA_encoders_chi1chi2/Encoders/', model))

All_states = {}
for state in os.listdir('/home/sarah/Downloads/AA_encoders/States'):
    All_states[state[-7:-4]] = np.loadtxt(os.path.join('/home/sarah/Downloads/AA_encoders_chi1chi2/States/', state))

for file in os.listdir('/home/sarah/Downloads/pdb_torsion_angles'):
    angles = pd.read_csv(os.path.join('/home/sarah/Downloads/pdb_torsion_angles/', file))
    con_seq = []
    drop = []
    tor_ls = []
    tor = angles.iloc[:,[5,7]]
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

    for i in range(angles.shape[0]):
        resi = angles.iloc[i,3]
        if resi not in ['UNK', 'GLY','ALA']:
            model = models[resi]
            states = All_states[resi]
            con_seq.append(concode[f'{resi}{discretize(tor_np[i,:], model, states)}'])
        else:
            con_seq.append(concode[resi])

    output_char = ''.join(con_seq)
    print(output_char)
    with open(f'/home/sarah/Downloads/pdb_con_seq_encoder_sep_chi1chi2/{file[7:14]}.txt','w') as f:
        f.write(output_char)
