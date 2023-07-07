import matplotlib.pyplot as plt
import pandas as pd 
import plotly.graph_objects as go
import os
from sklearn.cluster import DBSCAN 
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import numpy as np
import warnings
from plotly.express import scatter as px_scatter
from math import cos
from math import sin
warnings.filterwarnings('ignore')
os.chdir('/home/sarah/Downloads/Codes/')


colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'pink', 'brown', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'lavender', 'turquoise', 'darkgreen', 'tan', 'salmon', 'gold']


input_file = pd.read_csv('torsionlist.csv')
input_file = pd.DataFrame(input_file)
input_file_del = input_file[['resn', 'chi1', 'altchi1', 'chi2', 'altchi2', 'chi3', 'chi4', 'chi5']]
AAs = ['LYS']
for i in AAs:
    AA = input_file_del[input_file_del['resn'] == i]
    AA = AA.dropna(axis=1)
    data = AA.iloc[:,1:3]
    for a in range(0,2):
        anglecos = []
        for b in AA.iloc[:,a+1]:
            anglecos.append(cos(b/180 * 3.1415))
        data[f'{AA.columns[a+1]}cos'] = anglecos

        anglesin = []
        for b in AA.iloc[:,a+1]:
            anglesin.append(sin(b/180 * 3.1415))
        data[f'{AA.columns[a+1]}sin'] = anglesin
    print(data)
    data=data.to_numpy()
    print(data.shape)
    #db = DBSCAN(eps=10, min_samples=5).fit(data[:,2:])
    #labels = db.labels_
    #n_clusters = len(set(labels))-(1 if -1 in labels else 0)
    #cluster_id = KMeans(n_clusters).fit_predict(data[:,2:])
    cluster_id_sincos = KMeans(8).fit_predict(data[:,2:])

    #dbt = DBSCAN(eps=10, min_samples=5).fit(data[:,0:2])
    #labels = dbt.labels_
    #n_clusters = len(set(labels))-(1 if -1 in labels else 0)
    #cluster_idt = KMeans(n_clusters).fit_predict(data[:,0:2])
    cluster_id_chi12 = KMeans(8).fit_predict(data[:,0:2])
    fig1 = plt.scatter(x = data[:,0], y = data[:,1], c=cluster_id_sincos)
    plt.show()
    fig2 = plt.scatter(x = data[:,0], y = data[:,1], c=cluster_id_chi12)
    plt.show()