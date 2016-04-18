from __future__ import division
import os, argparse, random, math
from collections import Counter
from itertools import product
import pandas as pd

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'

class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print "Not in FASTA format."

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list

def gc(n):
    count = 0
    length = len(n)
    for i in n:
        if i == "C" or i == "G":
            count += 1
    GCrel = count / length
    return round(GCrel, 3)

def promotorToFeatures():
    fasta = mydir + 'DREAM6_ExPred_Promoters.fasta'
    class_test = classFASTA(fasta)
    OUT =  open(mydir+'DREAM6_ExPred_Promoters_Features.txt', 'w')
    # print sequence name, length, and GC-content to file
    for x in class_test.readFASTA():
        lenSeq = len(x[1])
        GCcontent = gc(x[1])
        print>> OUT, x[0], lenSeq, GCcontent
    OUT.close()

def mergeFeatureSeet():
    feature1 = pd.read_csv(mydir + 'DREAM6_ExPred_Promoters_Features.txt', \
        sep = ' ', header = None, names = ["Sequence", "Length", "GC"])
    feature2 = pd.read_csv(mydir + 'DREAM6_ExPred_PromoterActivities.txt', \
        sep = '\t', header = None, names = ["Sequence", "Expression"])
    featuresMerged = pd.merge(feature1, feature2, left_on = 'Sequence', right_on = 'Sequence')
    featuresMerged.to_csv(path_or_buf = mydir + 'DREAM6_ExPred_Promoters_Features_Activitis.txt', \
        sep = '\t')

mergeFeatureSeet()
