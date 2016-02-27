#Build GHMM to predict transmembrane domain of a given protein
Input:training set file containing 123 E.coli protein sequences and a test set file containg 2 transmembrane protein and 1 non-transmembrane protein
Process:PredictMO_ghmm.py
Output:features of protein sequences
##Build GHMM from training set


1) read and parse the protein sequences with features

2) generate length distribution for membrane domain,inside and outside domain

3) emission possibilty matrix for three hidden states

4) transition probability between the states
I: inside of cell
M: cell membrane (MI: from I; MO from o)
O: outside of cell 

5) initial probability of the states

6) implemented Viterbi algorithm to find the most likely hidden states for the observed sequences

##Command line
python PredictMO_ghmm.py -i ../data/tmptest.txt -o Predicted_Hidden_states.txt 
