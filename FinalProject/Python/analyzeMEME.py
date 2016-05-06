from __future__ import division
import os, argparse, random, math
from collections import Counter
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cross_validation import train_test_split, cross_val_score, LeaveOneOut,cross_val_predict
from sklearn.linear_model import LinearRegression, Lasso, ElasticNet
from sklearn.feature_selection import RFECV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn import cross_validation
import scipy.stats
from sklearn.preprocessing import Imputer

sns.set(style='whitegrid', context='notebook')

mydir = os.path.expanduser("~/github/I529-Group-Projects/FinalProject/")
IN = mydir + 'data/MEME_dataframe.txt'



def plotResiduals(y_train, y_test, y_train_pred, y_test_pred, filename = 'Fig1'):
    plt.scatter(y_train_pred, y_train_pred - y_train, \
        c = 'blue', marker = 'o', label = 'Training data')
    plt.scatter(y_test_pred, y_test_pred - y_test, \
        c = 'lightgreen', marker = 's', label = 'Test data')
    plt.xlabel('Predicted expression')
    plt.ylabel('Residuals')
    plt.legend(loc = 'upper left')
    plt.hlines(y = 0, xmin = -10, xmax = 50, lw =2, color = 'red')
    plt.xlim([-5, 5])
    plt.savefig(str(mydir) + 'figs/' +  filename + ".png")

def pairPlot(df, filename = 'Fig2'):
    cols = list(df.columns.values)
    sns_plot = sns.pairplot(df[cols], size=2.5)
    sns_plot.savefig(str(mydir) + 'figs/' + filename + ".png")

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')

def overallScore(obs, pred):
    pearsonsr = scipy.stats.pearsonr(pred, obs)
    CS = scipy.stats.chisquare(pred, f_exp = obs)
    spearmanr = scipy.stats.spearmanr(obs, pred)
    WSR = scipy.stats.wilcoxon(obs , y = pred)
    result =  (-1/4) *  np.log10(pearsonsr[1] * CS[1] * spearmanr[1] * WSR[1] )
    print "Pearsons r2 = " +str(pearsonsr[0])
    print "Pearsons p-value = " +str(pearsonsr[1])
    print "Spearman rho = " +str(spearmanr[0])
    print "Spearman p-value = " +str(spearmanr[1])
    print 'Chi-square test statistic = ' + str(CS[0])
    print 'Chi-square p-value = ' + str(CS[1])
    print 'Wilcox test statistic = ' + str(WSR[0])
    print 'Wilcox p-value = ' + str(WSR[1])
    print "overall score = " + str(result)
    return result

slr = LinearRegression()

df = pd.read_csv(IN, sep = '\t')
# Split predictor and response variables
X = df.drop(['Expression', 'Sequence'], axis=1)
X = X.ix[:, X.columns != 'Expression'].values
y = df['Expression'].values

imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
imp = imp.fit(X)
# Impute our data, then train
X_train_imp = imp.transform(X)

X_train_imp = pd.DataFrame(X_train_imp, \
    columns = ['Length(log10)', 'GC', '1_x','2_x','3_x','4_x','5_x','1_y','2_y',\
        '3_y','4_y','5_y','1','2','3','4','5'])

#X_train_imp = pd.DataFrame(X_train_imp, \
#    columns = ['Length(log10)', 'GC', '1_x','2_x','3_x','4_x','5_x','1_y','2_y',\
#        '3_y','4_y','5_y','1','2','3','4','5', '6_x', '6_y', '6'])

rfecv = RFECV(estimator=slr, step=1, cv=StratifiedKFold(y, 2))
rfecv = rfecv.fit(X_train_imp, y)
print("Optimal number of features : %d" % rfecv.n_features_)
print "RFECV ranking: " + str(rfecv.ranking_)
print "RFECV support: " + str(rfecv.support_)


'''Run with all features'''


print "Model 1"

predicted = cross_val_predict(slr, X_train_imp, y, cv=10)
fig, ax = plt.subplots()
ax.scatter(y, predicted)
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=4)
ax.set_xlabel('Measured gene expression')
ax.set_ylabel('Predicted gene expression')
plt.savefig(str(mydir) + 'figs/' +  'Fig3' + ".png")
plt.close()

overallScore(y, predicted)


'''Run with selected features.'''
X_train_imp = X_train_imp[['GC', '5_y']]
print "Model 2"
predicted1 = cross_val_predict(slr, X_train_imp, y, cv=10)
overallScore(y, predicted1)


# Figure 1

#pairPlot(X_train_imp)
# figure out how many features we need
# and what featuer it is




## plot cross validated observed predicted
# reshape X
#testX = np.reshape(X[:, 2:], (90, 6))





# pearsons

# spearman


#plt.figure()
#plt.xlabel("Number of features selected")
#plt.ylabel("Cross validation score (nb of correct classifications)")
#plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
#plt.savefig(str(mydir) + 'figs/' +  'RFECV' + ".png")




# Figure 2
X_train, X_test, y_train, y_test = train_test_split(X_train_imp, y, test_size = 0.3, \
    random_state = 0)
slr.fit(X_train, y_train)
y_train_pred = slr.predict(X_train)
y_test_pred = slr.predict(X_test)
plotResiduals(y_train, y_test, y_train_pred,  y_test_pred, filename = 'Fig2')
