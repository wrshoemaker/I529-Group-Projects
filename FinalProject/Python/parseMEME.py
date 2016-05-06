from __future__ import division
import pandas as pd
import numpy as np
import os, argparse, random, math
from collections import Counter
from itertools import product

mydir = os.path.expanduser("~/github/I529-Group-Projects/FinalProject/")


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


def getMotifCounts(IN = mydir + 'data/Meme_final_dna/meme.txt'):
    offset = 0
    if 'MW' in IN:
        offset = 1
    lengthDict = {}
    lengthCountDict = {}
    motifDict = {}
    pValueDict = {}
    for i, line in enumerate(open(IN)):
        line = line.strip().split()
        if len(line) == 0:
            continue

        if i in range(34, 79):
            for j in range(0, len(line), 3):
                lengthDict[line[j]] = int(line[j + 2])
                lengthCountDict[line[j]] = dict()
                pValueDict[line[j]] = dict()
                for k in range(1,6):
                    lengthCountDict[line[j]][str(k)] = []
                    pValueDict[line[j]][str(k)] = []

    for x in range(1, 6):
        motifDict[str(x)] = {}
        for key, value in lengthDict.iteritems():
            motifDict[str(x)][key] = float(0)

    motifNum = ''
    for i, line in enumerate(open(IN)):
        line = line.strip().split()
        if len(line) == 0 or line[0] == 'SEQUENCE' or line[-1] == '*****' or '(' in line[0] \
            or line[0] == 'Relative' or line[0] == 'Entropy' or line[0] == 'Stopped' \
            or line[0] == 'consensus' or line[0] == 'Sequence' or line[0] == 'Motif':
            continue
        if line[0] == 'MOTIF':
            motifNum = (line[1])
        if i > 150 and len(line) == 6 + offset:

            motifDict[motifNum][line[0]] += 1
            pValueDict[line[0]][motifNum].append(np.log10(float(line[2 + offset])))
            lengthCountDict[line[0]][motifNum].append(int(line[1 + offset]) / lengthDict[line[0]])
            #lengthDict[line[j]].append()
    for x in pValueDict:
        for y in pValueDict[x]:
            if len(pValueDict[x][y]) == 0:
                pValueDict[x][y] = np.nan
            else:
                pValueDict[x][y] = sum(pValueDict[x][y]) / len(pValueDict[x][y])
    for x in lengthCountDict:
        for y in lengthCountDict[x]:
            if len(lengthCountDict[x][y]) == 0:
                lengthCountDict[x][y] = np.nan
            else:
                lengthCountDict[x][y] = sum(lengthCountDict[x][y]) / len(lengthCountDict[x][y])
    return (motifDict, lengthCountDict, pValueDict)


mergeFeatureSet()

motifDictTest = getMotifCounts()
motifDictTestPandas =  pd.DataFrame.from_dict(motifDictTest[0])
motifDictTestPandas['Sequence'] = motifDictTestPandas.index

lengthCountDictTestPandas =  pd.DataFrame.from_dict(motifDictTest[1]).transpose()
lengthCountDictTestPandas['Sequence'] = lengthCountDictTestPandas.index
pValueDictTestPandas =  pd.DataFrame.from_dict(motifDictTest[2]).transpose()
pValueDictTestPandas['Sequence'] = pValueDictTestPandas.index

df = pd.read_csv(mydir + 'data/DREAM6_ExPred_Promoters_Features_Activitis.txt',\
    sep = '\t')
mergedDataFrames = pd.merge(df, motifDictTestPandas, on='Sequence', how='outer')
mergedDataFrames = pd.merge(mergedDataFrames, lengthCountDictTestPandas, on='Sequence', how='outer')
mergedDataFrames = pd.merge(mergedDataFrames, pValueDictTestPandas, on='Sequence', how='outer')


##### MW

motifDictTest_MW = getMotifCounts(IN = mydir + 'data/MEME_MW/meme.txt')
motifDictTestPandas_MW =  pd.DataFrame.from_dict(motifDictTest_MW[0])
motifDictTestPandas_MW['Sequence'] = motifDictTestPandas_MW.index

lengthCountDictTestPandas_MW =  pd.DataFrame.from_dict(motifDictTest_MW[1]).transpose()
lengthCountDictTestPandas_MW['Sequence'] = lengthCountDictTestPandas_MW.index
pValueDictTestPandas_MW =  pd.DataFrame.from_dict(motifDictTest_MW[2]).transpose()
pValueDictTestPandas_MW['Sequence'] = pValueDictTestPandas_MW.index

mergedDataFrames_MW = pd.merge(df, motifDictTestPandas_MW, on='Sequence', how='outer')
mergedDataFrames_MW = pd.merge(mergedDataFrames_MW, lengthCountDictTestPandas_MW, on='Sequence', how='outer')
mergedDataFrames_MW = pd.merge(mergedDataFrames_MW, pValueDictTestPandas_MW, on='Sequence', how='outer')
mergedDataFrames_MW.to_csv(path_or_buf = mydir + 'data/MEME_dataframe_MW.txt', \
    sep = '\t', index=False)

MW_to_merge =  mergedDataFrames_MW[['Sequence','1_x', '1_y', '1']]
MW_to_merge.columns = ['Sequence', '6_x', '6_y', '6']
#mergedDataFrames = pd.merge(mergedDataFrames, MW_to_merge, on='Sequence', how='outer')



mergedDataFrames.to_csv(path_or_buf = mydir + 'data/MEME_dataframe.txt', \
    sep = '\t', index=False)
