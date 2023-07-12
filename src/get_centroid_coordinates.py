import pandas as pd 
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from math import cos,sin
import os
import warnings

warnings.filterwarnings('ignore')
os.chdir('/home/sarah/Downloads/Codes/')


coordinates = dict()
input_file = pd.read_csv('all_torsion_angles.csv')
input_file = pd.DataFrame(input_file)
input_file_del = input_file[['resn', 'chi1', 'chi2']]
AAs = ['ARG','ASN','ASP','CYS','GLU', 'GLN','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

for x,i in enumerate(AAs):
    print(x)
    AA = input_file_del[input_file_del['resn'] == i]
    AA = AA[[i for i in AA.columns if AA[[i]].isnull().all()[0]==False]]
    AA = AA.dropna(axis=0)
    data = AA.iloc[:,1:]
    for a in range(len(AA.columns)-1):
        anglecos = []
        for b in AA.iloc[:,a+1]:
            anglecos.append(cos(b/180 * 3.1415))
        data[f'{AA.columns[a+1]}cos'] = anglecos

        anglesin = []
        for b in AA.iloc[:,a+1]:
            anglesin.append(sin(b/180 * 3.1415))
        data[f'{AA.columns[a+1]}sin'] = anglesin
    data=data.to_numpy()
    cluster_id_sincos = KMeans(3).fit(data[:,len(AA.columns)-1:len(AA.columns)+3])
    cluster_coordinates = cluster_id_sincos.cluster_centers_

    coordinates[f'centroids{i}'] = cluster_coordinates


print(coordinates)
np.save('cluster_centroids.npy', coordinates)