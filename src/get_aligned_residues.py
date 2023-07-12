import re
import numpy as np
import os
import pandas as pd
from tqdm import tqdm

os.chdir('/home/sarah/Downloads/Codes')

def parse_cigar(cigar_string):
    ref, query = 0, 0
    matches = []

    for cnt, action in re.findall('([0-9]*)([IDMP])', cigar_string):
        cnt = int(cnt)

        if action == 'D':
            ref += cnt
        elif action == 'I':
            query += cnt
        elif action == 'M':  # use only Perfect matches
            ref += cnt
            query += cnt
        elif action == 'P':
            matches += [(ref + i, query + i) for i in range(cnt)]
            ref += cnt
            query += cnt
        else:
            raise ValueError(f'Action {action}')
    return np.array(matches)


tor_angles = {}
for a, name in enumerate(os.listdir('/home/sarah/Downloads/pdb_torsion_angles')):
    print(a)
    with open(os.path.join('/home/sarah/Downloads/pdb_torsion_angles/', name)) as file:
        file = pd.read_csv(file)
        tor_angles[name[7:-4]] = file.iloc[:,[3,5,6,7,8,9,10,11]]


pair_file = open('alignment_cut.txt')

alignment_list = list()
for line in tqdm(pair_file):
    sid1, sid2, cigar_string = line.rstrip('\n').split()
    idx_1, idx_2 = parse_cigar(cigar_string).T
    sid1 = tor_angles[sid1].to_numpy()
    sid2 = tor_angles[sid2].to_numpy()
    alignment_list.append(np.hstack((sid1[idx_1], sid2[idx_2])))
alignment_array = np.vstack(alignment_list)
alignment_df = pd.DataFrame(alignment_array, columns=['resi1','chi1_1','altchi1_1','chi2_1','altchi2_1','chi3_1','chi4_1','chi5_1','resi2','chi1_2','altchi1_2','chi2_2','altchi2_2','chi3_2','chi4_2','chi5_2'])
alignment_df.to_csv('alignment_cut_list.csv')



    

