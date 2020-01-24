# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:36:36 2020

@author: miseminger
"""

import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import ExpressionModel


x = np.array([0, 1.25, 2.5, 5, 10, 20]) #protein concentration (nM)
y = np.array([0.1475, 0.334, 0.444, 0.532, 0.598, 0.6615]) #mean absorbance
y_err = np.array([0.0175, 0.003, 0.013, 0.003, 0.013, 0.0235]) #standard error of the mean


fit_model = ExpressionModel('B * x**h / (Kd**h + x**h)', independent_vars=['x']) #Hill function
params = fit_model.make_params(B=0.7, h=0.8, Kd=1.6) #include guesses for parameters
params['B'].min = 0.5 #set minimum values for parameters
params['h'].min = 0.5
params['Kd'].min = 1


result = fit_model.fit(y, params, x=x, weights=1.0/y_err)
print(result.fit_report())


#plot

x_mod = np.linspace(0,25,num=1000)

B_best = result.params['B'].value
h_best = result.params['h'].value
Kd_best = result.params['Kd'].value
y_best = B_best * x_mod**h_best / (Kd_best**h_best + x_mod**h_best)

B_init = params['B'].value
h_init = params['h'].value
Kd_init = params['Kd'].value
y_init = B_init * x_mod**h_init / (Kd_init**h_init + x_mod**h_init)

plt.errorbar(x, y, yerr=y_err, capsize=3, fmt = '.', color='b')
plt.plot(x_mod, y_init, 'k--', label='initial fit')
plt.plot(x_mod, y_best, 'r-', label='best fit')
plt.legend(loc='best')
plt.xlabel('protein concentration (nM)')
plt.ylabel('OD 450')
plt.show()
