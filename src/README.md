angle_distributions.py:
Angle distributions (histograms)


get_torsion.py: atomic coordinates to torsion angles 
from	https://gist.github.com/lennax/0f5f65ddbfa278713f58


get_centroid_coordinates.py: 
input: a csv file containing all residues and their angles
ouput: a npy file containing coordinates of three clusters for each AA except GLY and ALA


alphabet_conversion.py:
loads the centroid coordinates
creates an alphabet
loads every sequence individually and creates a file with the converted sequence.


combine_con_seqs.py:
input: converted sequences in a directory
output: a csv file containing ids and sequences (for substitution matrix reconstruction)


create_submat.py:
substitution matrix reconstruction
from	https://github.com/steineggerlab/foldseek-analysis/blob/main/training/create_submat.py


get_aligned_residues.py:
prepares a csv file of aligned pairs of residues, required for training vqvae.


prepare_alignment_vectors.py:
input: csv file of residues and their angles
output: vectors for training the autoencoder 
(first 20 elements specify the AA, the rest are sin/cos of the angles of interest, NA values:2)


train_vqvae.py:
line 131 should be changed based on the angles that are used.
from	https://github.com/steineggerlab/foldseek-analysis/blob/main/training/train_vqvae.py


vqvae_convert_seqs.py:
loads encoder and decoder
converts sequences and create a separate file for each of them


prepare_alignment_separate.py:
prepares a separate set of alignment vectors for each AA


train_vqvae_separate.py:
creates a separate model for each amino acid


vqvae_convert_seqs_separate.py:
converts the sequences based on separate models for each amino acid




