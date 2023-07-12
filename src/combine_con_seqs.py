import os
import pandas as pd
import csv

data = pd.DataFrame()
l_sid = []
l_sequence = []
for f in os.listdir('/home/sarah/Downloads/pdb_con_seq_encoder_sep_chi1chi2'):
    l_sid.append(f[:-4])
    print(f[:-4])
    with open(os.path.join('/home/sarah/Downloads/pdb_con_seq_encoder_sep_chi1chi2/', f)) as seq:
        l_sequence.append(seq.read())
data['sid'] = l_sid
data['sequence'] = l_sequence


os.chdir('/home/sarah/Downloads/Codes')
data.to_csv('combined_con_seq_encoder_sep_chi1chi2.csv', lineterminator='\n')