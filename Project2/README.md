#Build GHMM to predict transmembrane domain of a given protein
Input:training set file containing 123 E.coli protein sequences and a test set file containg 2 transmembrane protein and 1 non-transmembrane protein
Process:PredictMO_ghmm.py
Output:features of protein sequences
##Build GHMM from training set


1)read and parse the protein sequences with features
i: inside of cell
m: cell membrane
o: outside of cell
)generate length distribution for membrane domain,inside and outside domain

3)emission possibilty matrix for three hidden states

#predict transmembrane domain of given protein sequence in fasta format
