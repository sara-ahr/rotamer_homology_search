import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


file = pd.read_csv('/home/sarah/Downloads/Codes/all_torsion_angles.csv', low_memory=False)
torsion_angles = pd.DataFrame(file)
print(torsion_angles.iloc[2,3])

font = {'family': 'sans',
        'color':  'black',
        'weight': 'normal',
        'size': 14,
        }

degree_sign = u'\N{DEGREE SIGN}'

amino_acids = {
    "ALA": [],  # Alanine
    "ARG": ["chi1", "chi2", "chi3", "chi4", "chi5"],  # Arginine
    "ASN": ["chi1", "chi2"],  # Asparagine
    "ASP": ["chi1", "chi2", "altchi2"],  # Aspartic Acid
    "CYS": ["chi1"],  # Cysteine
    "GLN": ["chi1", "chi2", "chi3"],  # Glutamine
    "GLU": ["chi1", "chi2", "chi3"],  # Glutamic Acid
    "GLY": [],  # Glycine
    "HIS": ["chi1", "chi2"],  # Histidine
    "ILE": ["chi1", "chi2"],  # Isoleucine
    "LEU": ["chi1", "chi2", "altchi2"],  # Leucine
    "LYS": ["chi1", "chi2", "chi3", "chi4"],  # Lysine
    "MET": ["chi1", "chi2", "chi3"],  # Methionine
    "PHE": ["chi1", "chi2", "altchi2"],  # Phenylalanine
    "PRO": ["chi1", "chi2"],  # Proline
    "SER": ["chi1"],  # Serine
    "THR": ["chi1"],  # Threonine
    "TRP": ["chi1", "chi2"],  # Tryptophan
    "TYR": ["chi1", "chi2", "altchi2"],  # Tyrosine
    "VAL": ["chi1", "altchi1"],  # Valine
}


def get_dist(AA, angle, torsion_angles):
    angle_list = []
    angles = list(torsion_angles.columns)
    angle_idx = angles.index(angle)
    for a,i in enumerate(torsion_angles['resn']):
        if i == AA:
            angle_list.append(torsion_angles.iloc[a, angle_idx])
    fig = plt.hist(angle_list, bins=range(-180, 180, 10), color= 'darkcyan')
    plt.ylabel('Frequency')
    plt.xlabel(f'Angle ({degree_sign})')
    plt.xticks(np.arange(-180,180,45))
    plt.title(f'Distribution of {angle} in {AA} residues') 
    return(fig)

for i in amino_acids.keys():
    for j in torsion_angles.columns[-7:]:
        if j in amino_acids[i]:
            fig = get_dist(i,j,torsion_angles)
            plt.savefig(f'/home/sarah/Documents/Report/Angle_distributions/{i}{j}')

        
