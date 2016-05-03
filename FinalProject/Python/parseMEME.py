from __future__ import division
import pandas as pd
import numpy as np
import os, argparse, random, math
from collections import Counter
from itertools import product

mydir = os.path.expanduser("~/github/I529-Group-Projects/FinalProject/")
IN = mydir + 'data/Meme_final/meme.txt'

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
                    print "Inproper file format."
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
    fasta = mydir + 'data/DREAM6_ExPred_Promoters.fasta'
    class_test = classFASTA(fasta)
    OUT =  open(mydir+'data/DREAM6_ExPred_Promoters_Features.txt', 'w')
    for x in class_test.readFASTA():
        lenSeq = np.log10(len(x[1]))
        GCcontent = gc(x[1])
        print>> OUT, x[0], lenSeq, GCcontent
    OUT.close()

def mergeFeatureSet():
    feature1 = pd.read_csv(mydir + 'data/DREAM6_ExPred_Promoters_Features.txt', \
        sep = ' ', header = None, names = ["Sequence", "Length(log10)", "GC"])
    feature2 = pd.read_csv(mydir + 'data/DREAM6_ExPred_PromoterActivities.txt', \
        sep = '\t', header = None, names = ["Sequence", "Expression"])
    featuresMerged = pd.merge(feature1, feature2, left_on = 'Sequence', right_on = 'Sequence')
    featuresMerged.to_csv(path_or_buf = mydir + 'data/DREAM6_ExPred_Promoters_Features_Activitis.txt', \
        sep = '\t', index=False)


def getMotifCounts():
    lengthDict = {}
    motifDict = {}
    for i, line in enumerate(open(IN)):
        line = line.strip().split()
        if len(line) == 0:
            continue

        if i in range(34, 79):
            for j in range(0, len(line), 3):
                lengthDict[line[j]] = int(line[j + 2])
    for x in range(1, 6):
        motifDict[str(x)] = {}
        for key, value in lengthDict.iteritems():
            motifDict[str(x)][key] = 0

    motifNum = ''
    for i, line in enumerate(open(IN)):
        line = line.strip().split()
        if len(line) == 0 or line[0] == 'SEQUENCE' or line[-1] == '*****':
            continue
        if line[0] == 'Motif':
            motifNum = (line[1])
        if i > 150 and len(line) == 6:
            motifDict[motifNum][line[0]] += 1
    return motifDict

mergeFeatureSet()
motifDictTest = getMotifCounts()
motifDictTestPandas =  pd.DataFrame.from_dict(motifDictTest)
motifDictTestPandas['Sequence'] = motifDictTestPandas.index


df = pd.read_csv(mydir + 'data/DREAM6_ExPred_Promoters_Features_Activitis.txt',\
    sep = '\t')
mergedDataFrames = pd.merge(df, motifDictTestPandas, on='Sequence', how='outer')

mergedDataFrames.to_csv(path_or_buf = mydir + 'data/MEME_dataframe.txt', \
    sep = '\t', index=False)
