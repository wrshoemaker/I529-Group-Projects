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
from sklearn.metrics import mean_squared_error, r2_score
from sklearn import cross_validation
sns.set(style='whitegrid', context='notebook')

mydir = os.path.expanduser("~/github/I529-Group-Projects/FinalProject/")
IN = mydir + 'data/Meme_final/meme.txt'



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
    plt.savefig(str(mydir) + 'figs/' +  filename + ".png")

def pairPlot(df):
    cols = ['Length(log10)','GC','Expression', '1', '2', '3', '4', '5']
    sns_plot = sns.pairplot(df[cols], size=2.5)
    sns_plot.savefig(str(mydir) + 'figs/' + "output.png")



df = pd.read_csv(mydir + 'data/MEME_dataframe.txt',\
    sep = '\t')

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


# cross validated


# cross_val_predict returns an array of the same size as `y` where each entry
# is a prediction obtained by cross validated:
predicted = cross_val_predict(slr, X, y, cv=20)
fig, ax = plt.subplots()
ax.scatter(y, predicted)
ax.plot([y.min(), y.max()], [y.min(), y.max()], 'k--', lw=4)
ax.set_xlabel('Measured')
ax.set_ylabel('Predicted')
plt.savefig(str(mydir) + 'figs/' +  'testCV' + ".png")


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

plotResiduals(y_train, y_test, y_train_pred,  y_test_pred, filename = 'residuals1')
