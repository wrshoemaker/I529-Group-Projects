from __future__ import division
import pandas as pd
import os

mydir = os.path.expanduser("~/github/I529-Group-Projects/FinalProject/")
IN = mydir + 'data/Meme_final/meme.txt'


# motif dict structure
# motif# { sequenceName {count, (start position(s))  } }

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
motifDictTest = getMotifCounts()
motifDictTestPandas =  pd.DataFrame.from_dict(motifDictTest)
motifDictTestPandas['Sequence'] = motifDictTestPandas.index


df = pd.read_csv(mydir + 'data/DREAM6_ExPred_Promoters_Features_Activitis.txt',\
    sep = '\t')
mergedDataFrames = pd.merge(df, motifDictTestPandas, on='Sequence', how='outer')

mergedDataFrames.to_csv(path_or_buf = mydir + 'data/MEME_dataframe.txt', \
    sep = '\t', index=False)
