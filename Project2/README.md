#Build GHMM to predict transmembrane domain of a given protein
We have been provdied a training set file containing 123 E.coli protein sequences and a test set file containg 2 transmembrane protein and 1 non-transmembrane protein

##Build GHMM from training set
1)read and parse the protein sequences with features
i: inside of cell
m: cell membrane
o: outside of cell
2)generate length distribution for membrane domain and non-membrane domain

3)emission possibilty matrix for two hidden states

###predict transmembrane domain of given protein sequence in fasta format
