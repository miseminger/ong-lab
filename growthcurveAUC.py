# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:07:37 2020

@author: miseminger
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import trapz
from scipy.integrate import simps
from sklearn.metrics import auc
import pandas as pd

#sample data that is missing weekends
x = np.append(np.arange(0,5), np.arange(7,12))
viable_cells = np.array([201069.9775, 384066.7396, 893906.1458, 2618218.583, 7317207.833, 15432305, 13953402.33, 10717670.33, 5858551.667, 1058330.896])
y_err = np.array([0, 13906.56944, 43355.77778, 241728.8889, 477322.5556, 1519088.667, 1295765.111, 655244.8889, 1774314.778, 580394.8056])



#combine x and y into a pandas dataframe with x as the index
df = pd.DataFrame(data=viable_cells, index=x, columns=['viable_cells'])
#reindex it to fill in empty numbers in x
newindex = np.arange(0,viable_cells.shape[0]+2)
df = df.reindex(newindex)
#interpolate missing data points so you can estimate the AUC
df_interp = df.interpolate(how='linear')

#plot interpolated data against initial data 
plt.errorbar(x, viable_cells, yerr=y_err, fmt="-o")
plt.scatter(np.array(df_interp.index), df_interp['viable_cells'], color='r')

#get area under the curve
area = trapz(df_interp['viable_cells'], dx=1) 
print("numpy.trapz area =", area)

area = simps(df_interp['viable_cells'], dx=1) #can also work with unevenly spaced data: can feed an x and y array
print("scipy.integrate.simps area =", area)

area = auc(np.array(df_interp.index), df_interp['viable_cells'])  #auc actually assumes monotonicity, which this data doesn't have
print("sklearn.metrics.auc area =", area)

'''
Output:

numpy.trapz area = 80554541.89815
scipy.integrate.simps area = 80969810.02590832
sklearn.metrics.auc area = 80554541.89815
'''
