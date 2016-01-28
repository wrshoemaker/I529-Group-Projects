import pandas as pd
from scipy.stats import gaussian_kde
import numpy as np
import  matplotlib.pyplot as plt


test = pd.read_csv("./test.txt", header=True, delim_whitespace=True)

def get_kdens_choose_kernel(xlist,expand, kernel = 0.5):
    """ Finds the kernel density function across a vector of values """
    xlist = xlist[np.logical_not(np.isnan(xlist))]
    density = gaussian_kde(xlist)
    n = len(xlist)
    if expand == False:
        xs = np.linspace(min(xlist),max(xlist),n)
    else:
        xs = np.linspace(min(xlist - expand),max(xlist + expand),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

LL = test[test.columns[2]].values

kde_data = get_kdens_choose_kernel(LL, 0.4, kernel = 0.4)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(kde_data[0], kde_data[1])
plt.xticks(fontsize = 6) # work on current fig
plt.yticks(fontsize = 6)
plt.ylabel('Probability Density', fontsize = 14)
plt.xlabel('Likelihood', fontsize = 14)

plt.tight_layout()
plt.savefig("test_plot.png", bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.xscale()
plt.close()
