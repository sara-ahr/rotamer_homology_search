import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
from sklearn.cluster import KMeans
from math import cos, sin
import warnings
import os

warnings.filterwarnings('ignore')
os.chdir('/home/sarah/Downloads/Codes/')


coordinates = np.load('cluster_centroids.npy', allow_pickle = True).item()

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
    concode['miss_info'] = chr(AAnum)
else:    
    while len(chr(AAnum)) > 1:
        AAnum += 1
    concode['miss_info'] = chr(AAnum)

for a in ['GLY','ALA']:
    AAnum += 1
    if len(chr(AAnum)) == 1:
        concode[a] = chr(AAnum)
    else:
        while len(chr(AAnum)) > 1:
            AAnum += AAnum
        concode[a] = chr(AAnum)
print(concode)


    
for f in os.listdir('/home/sarah/Downloads/pdb_torsion_angles'):
    inputname = os.path.join('/home/sarah/Downloads/pdb_torsion_angles/', f)
    input_file = pd.read_csv(inputname)
    input_file = pd.DataFrame(input_file)
    input_file_del = input_file[['resn', 'chi1', 'chi2']]
    
    sequence =[]
    for i in range(len(input_file_del)):
        if input_file_del.iloc[i,0] not in ['GLY','ALA']:
            query = np.array(input_file_del.iloc[i,1:])
            query = query[~pd.isnull(query)]
            query = query[:2]/180 * 3.1415
            if len(query)>0: 
                try:
                    query = np.array([cos(query[0]), sin(query[0]), cos(query[1]), sin(query[1])])
                except IndexError:
                    query = np.array([cos(query[0]), sin(query[0])])
                try:
                    label = np.argmin(np.linalg.norm(coordinates[f'centroids{input_file_del.iloc[i,0]}']-np.array([query, query, query]), axis=1))
                    sequence.append(concode[f'{input_file_del.iloc[i,0]}{label}'])
                except ValueError:
                    sequence.append(concode['miss_info'])
            else:    
                sequence.append(concode['miss_info'])
        else:
            sequence.append(concode[input_file_del.iloc[i,0]])
        

    name = inputname.split('/')[-1]
    name = name.split('.')[:-1]
    name = '.'.join(name)
    name = str(name[7:])
    os.chdir('/home/sarah/Downloads/pdb_con_seq_KMeans_final/')
    with open(f'con_seq_{name}', 'w') as f:
        f.write(''.join(sequence))
    print(f'con_seq_{name}','completed')
