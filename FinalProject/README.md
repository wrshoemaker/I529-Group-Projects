# I529: Final Project

### Authors 

Will Shoemaker, Xiangyu Yao and Meng Wu

### Run the analysise 

To generate the output from MEME, install MEME and run the following command:

	meme DREAM6_ExPred_Promoters.fasta -o Meme_final_dna -dna -mod anr -nmotifs 5 -minw 4 -maxw 15

There are two scripts to complete the analysis. 

First, run the following to generate the data matrix:

	python parseMEME.py

Second, run the following to get the statistical output and the figures:

	python analyzeMEME.py

This completes the analysis.
	
