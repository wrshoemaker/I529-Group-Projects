from __future__ import division
import os, argparse, random, math
from collections import Counter
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import LinearRegression, Lasso, ElasticNet
from sklearn.metrics import mean_squared_error, r2_score
from sklearn import cross_validation
sns.set(style='whitegrid', context='notebook')

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
        lenSeq = np.log10(len(x[1]))
        GCcontent = gc(x[1])
        print>> OUT, x[0], lenSeq, GCcontent
    OUT.close()

def mergeFeatureSet():
    feature1 = pd.read_csv(mydir + 'DREAM6_ExPred_Promoters_Features.txt', \
        sep = ' ', header = None, names = ["Sequence", "Length(log10)", "GC"])
    feature2 = pd.read_csv(mydir + 'DREAM6_ExPred_PromoterActivities.txt', \
        sep = '\t', header = None, names = ["Sequence", "Expression"])
    featuresMerged = pd.merge(feature1, feature2, left_on = 'Sequence', right_on = 'Sequence')
    featuresMerged.to_csv(path_or_buf = mydir + 'DREAM6_ExPred_Promoters_Features_Activitis.txt', \
        sep = '\t', index=False)

def plotResiduals(y_train, y_test, y_train_pred, y_test_pred, filename = 'residuals'):
    plt.scatter(y_train_pred, y_train_pred - y_train, \
        c = 'blue', marker = 'o', label = 'Training data')
    plt.scatter(y_test_pred, y_test_pred - y_test, \
        c = 'lightgreen', marker = 's', label = 'Test data')
    plt.xlabel('Predicted values')
    plt.ylabel('Residuals')
    plt.legend(loc = 'upper left')
    plt.hlines(y = 0, xmin = -10, xmax = 50, lw =2, color = 'red')
    plt.xlim([-5, 5])
    plt.savefig(str(mydir[:-6]) + '/figs/' +  filename + "png")

df = pd.read_csv(mydir + 'DREAM6_ExPred_Promoters_Features_Activitis.txt',\
    sep = '\t')
#cols = ['Length(log10)','GC','Expression']
#sns_plot = sns.pairplot(IN[cols], size=2.5)
#sns_plot.savefig(str(mydir[:-6]) + '/figs/' + "output.png")

X = df.iloc[:, 1:-1].values
y = df['Expression'].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, \
    random_state = 0)
slr = LinearRegression()
slr.fit(X_train, y_train)
y_train_pred = slr.predict(X_train)
y_test_pred = slr.predict(X_test)

print('MSE train: %.3f, test: %.3f' % (
        mean_squared_error(y_train, y_train_pred),
        mean_squared_error(y_test, y_test_pred)))
print('R^2 train: %.3f, test: %.3f' % (
        r2_score(y_train, y_train_pred),
        r2_score(y_test, y_test_pred)) )

# Lets try looking at the data with LASSO and Elastic Net
EN = ElasticNet(alpha = 1.0, l1_ratio = 0.5)
lasso = Lasso(alpha = 1.0)

scores1 = cross_validation.cross_val_score( \
    EN, X_train, y_train, cv=20, scoring='r2')
# This will print the mean of the list of errors that were output and
print "Elastic Net"
print scores1.mean()
print scores1.std()

scores2 = cross_validation.cross_val_score( \
    lasso, X_train, y_train, cv=20, scoring='r2')

print "Lasso"
print scores2.mean()
print scores2.std()
