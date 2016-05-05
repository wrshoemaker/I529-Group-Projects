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
sns.set(style='whitegrid', context='notebook')

mydir = os.path.expanduser("~/github/I529-Group-Projects/FinalProject/")
IN = mydir + 'data/Meme_final/meme.txt'



def plotResiduals(y_train, y_test, y_train_pred, y_test_pred, filename = 'Fig1'):
    plt.scatter(y_train_pred, y_train_pred - y_train, \
        c = 'blue', marker = 'o', label = 'Training data')
    plt.scatter(y_test_pred, y_test_pred - y_test, \
        c = 'lightgreen', marker = 's', label = 'Test data')
    plt.xlabel('Predicted values')
    plt.ylabel('Residuals')
    plt.legend(loc = 'upper left')
    plt.hlines(y = 0, xmin = -10, xmax = 50, lw =2, color = 'red')
    plt.xlim([-5, 5])
    plt.savefig(str(mydir) + 'figs/' +  filename + ".png")

def pairPlot(df, filename = 'Fig2'):
    cols = ['Length(log10)','GC','Expression', '1', '2', '3', '4', '5']
    sns_plot = sns.pairplot(df[cols], size=2.5)
    sns_plot.savefig(str(mydir) + 'figs/' + filename + ".png")


df = pd.read_csv(mydir + 'data/MEME_dataframe.txt',\
    sep = '\t')
# Split predictor and response variables
X = df.drop(['Expression', 'Sequence'], axis=1)
X = X.ix[:, X.columns != 'Expression'].values
y = df['Expression'].values
EN = ElasticNet(alpha = 1.0, l1_ratio = 0.5)
slr = LinearRegression()
# Figure 1
#pairPlot(df)

# figure out how many features we need
# and what featuer it is

rfecv = RFECV(estimator=slr, step=1, cv=StratifiedKFold(y, 2))
rfecv = rfecv.fit(X, y)
print("Optimal number of features : %d" % rfecv.n_features_)
print "RFECV ranking: " + str(rfecv.ranking_)
print "RFECV support: " + str(rfecv.support_)



## plot cross validated observed predicted
# reshape X
#testX = np.reshape(X[:, 2:], (90, 6))
predicted = cross_val_predict(slr, X, y, cv=20)
fig, ax = plt.subplots()
ax.scatter(y, predicted)
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=4)
ax.set_xlabel('Measured gene expression')
ax.set_ylabel('Predicted gene expression')
plt.savefig(str(mydir) + 'figs/' +  'Fig3' + ".png")


print scipy.stats.pearsonr(y, predicted)
print scipy.stats.chisquare(predicted, f_exp = y)
print scipy.stats.spearmanr(y, predicted)

# pearsons

# spearman


#plt.figure()
#plt.xlabel("Number of features selected")
#plt.ylabel("Cross validation score (nb of correct classifications)")
#plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
#plt.savefig(str(mydir) + 'figs/' +  'RFECV' + ".png")




# Figure 2
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, \
    random_state = 0)
plotResiduals(y_train, y_test, y_train_pred,  y_test_pred, filename = 'Fig2')







scores = cross_val_score(EN, X, y, cv=20, scoring='r2')
print scores.mean() *-1
print scores.std()

scores1 = cross_val_score(EN, X[:,3:4], y, cv=20, scoring='r2')
print scores1.mean() *-1
print scores1.std()













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


# cross validated recursive feature plot


# Create the RFE object and compute a cross-validated score.







#selector = RFECV(slr, step=1, cv=StratifiedKFold(y, 2))
#selector = selector.fit(X, y)
#print selector.support_
#print selector.ranking_

# cross_val_predict returns an array of the same size as `y` where each entry
# is a prediction obtained by cross validated:


# Lets try looking at the data with LASSO and Elastic Net
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

plotResiduals(y_train, y_test, y_train_pred,  y_test_pred, filename = 'residuals1')
