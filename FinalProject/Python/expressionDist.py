from __future__ import division
import pylab as plt
import pandas as pd
import numpy as np
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde


def CV_KDE(oneD_array, expand = 100):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 1.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    # add nothing to the end of grid and pdf so you can get a nice looking kde
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    #x_grid = np.linspace(min(x_grid - expand),max(x_grid + expand),100)
    #pdf = np.lib.pad(pdf, (100,100), 'constant', constant_values=(0, 0))
    return_tuple = (x_grid, pdf, kde.bandwidth)
    print kde.bandwidth
    return return_tuple

def get_kdens_choose_kernel(xlist,expand, kernel = 2):
    """ Finds the kernel density function across a vector of values """
    xlist = xlist[np.logical_not(np.isnan(xlist))]
    density = gaussian_kde(xlist)
    n = len(xlist)
    #n = 10
    if expand == False:
        xs = np.linspace(min(xlist),max(xlist),n)
    else:
        xs = np.linspace(min(xlist - expand),max(xlist + expand),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

# import data
#fig, ax = plt.subplots()

importData = pd.read_csv('../data/DREAM6_ExPred_PromoterActivities.txt', sep = '\t', header = None)
expression = importData[importData.columns[1]].values
#returnKDE = CV_KDE(expression)
#plt.hist(expression, 30, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
returnKDE = get_kdens_choose_kernel(expression, False, kernel = 0.25 )


plt.plot(returnKDE[0], returnKDE[1],color = 'b', linestyle = '-', label="N = 1000, B = 1")
#plt.fill_between(returnKDE[0], 0,returnKDE[1])
plt.ylabel('Probability', fontsize = 18, rotation = 90)
plt.xlabel('Expression level', fontsize = 18)
output = '../figs/testKDE.png'
plt.savefig(output, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.xscale()
plt.close()


plt.hist(expression, 30, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
output = '../figs/testHist.png'
plt.savefig(output, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
