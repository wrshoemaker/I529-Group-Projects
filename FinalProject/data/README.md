#generate training set and test set
fasta-subsample <input> 81 -rest test.fasta > training.fasta

#discover motifs
meme training.fasta -o Meme`_`out -dna 
