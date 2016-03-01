#Implement Felsenstein's pruning algorithm on three aligned sequence

##Input
Human		AGCTTC</br>
Chimpanzee	AGTTGC</br>
Gorilla		ACTTGC</br>

Phylogenetic tree for these three sequences in Newick Standard format.</br>
((Human:0.3,Chimpanzee:0.2):0.1,Gorilla:0.3)</br>

##Output
The probability to observe these sequences 

## Scheme of program
### Calculate the substitution matrix 
Given the base substitution matrix P(0.1), we need to compute the other matrix P(0.2) and P(0.3).

### Pruning algorithm to calcualte probabilty for each column in the alignment
Recurring calculation of probability of subtree rooted in Xi.

### Forward algorithm to sum probabilities of each column
