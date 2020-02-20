# -*- coding: utf-8 -*-

"""

Created on Sat Jan 18 19:07:03 2020

@author: madeline

"""


import pandas as pd

import matplotlib.pyplot as plt

import numpy as np

from lmfit.models import LinearModel

from datetime import datetime

from scipy.integrate import simps



#fill out these data things before you plot!

filename = r'E:\\cedexjan1to29_addtoflip.xlsx' #name of Excel file to read from

sheet = 'p14.4' #sheet name if applicable (0 is default)

start_row = 6 - 2

end_row = 12 - 2

linearpoints = 3 #number of points in linear phase to fit

first_linearpoint = 0 #first point is 0

title = 'P14.4 Jan 21' #plot title





#read in the data from Excel!

skip_rows = start_row - 1 #number of rows to skip

df = pd.read_excel(filename, sheet_name=sheet)



num_data_points = end_row - start_row

days = df['Date finished'][start_row:end_row+1]

viability = df['Viability'][start_row:end_row+1]

#viability_err = df['Mean_Viability'][start_row:end_row+1]
viability_err = 0

viable = df['Viable Cell Conc.'][start_row:end_row+1]/1000000

#viable_err = df['Mean_Viable'][start_row:end_row+1]/1000000
viable_err = 0

#translate dates to hour differences

hours = np.zeros((len(days)))

i=0

for j in np.arange(start_row, end_row+1):

    first_day = datetime.strptime(days[start_row], '%m/%d/%Y %I:%M:%S %p')

    new_day = datetime.strptime(days[j], '%m/%d/%Y %I:%M:%S %p')

    diff = new_day - first_day

    diffhours = int(diff.total_seconds()//(3600))

    hours[i] = diffhours/24 #this actually makes into days (floating point integers)

    if i != len(days) - 1:

        i=i+1

    

#set basis for figure names

fignamebase = '/Users/miseminger/Desktop/plots/' + title + str(start_row) + '_' + str(end_row)





#plot % viability

fig = plt.figure(1)

plt.errorbar(hours,viability, yerr=viability_err, color='g', fmt='-o')       

plt.xlabel('days')

plt.ylabel('% viable')

plt.ylim(0,110)

figname = fignamebase + '_percentviable.png'

plt.savefig(figname, dpi=199)





#plot viable cell concentration in millions/ml

fig = plt.figure(2)

plt.errorbar(hours,viable, yerr=viable_err, fmt='-o')

plt.xlabel('days')

plt.ylabel('viable cells (millions/ml)')

figname = fignamebase + '_viablecellconc.png'

plt.savefig(figname, dpi=199)





#plot semilog of viable cell concentration

fig = plt.figure(3)



viable_log = np.log(np.array(viable).astype(float))



mod = LinearModel()

pars = mod.guess(viable_log[first_linearpoint:first_linearpoint + linearpoints], x=hours[first_linearpoint:first_linearpoint + linearpoints])

init = mod.eval(pars, x=hours[first_linearpoint:first_linearpoint + linearpoints])

out = mod.fit(viable_log[first_linearpoint:first_linearpoint+linearpoints], pars, x=hours[first_linearpoint:first_linearpoint + linearpoints])



slope = out.params['slope'].value

intercept = out.params['intercept'].value

y = slope * hours + intercept    

doubling_time = int(np.log(2)/slope*24) #doubling time in hours



plt.plot(hours, y, 'r-', label='linear fit')

plt.errorbar(hours, viable_log, yerr=viable_err/viable, label='viable cells', fmt='-o')  

plt.xlabel('days')

plt.ylabel('ln(cells/ml)')

plt.legend(loc='upper left')

plt.ylim(np.min(viable_log)-1, np.max(viable_log)+1)

figname = fignamebase + '_semilogviable_' + str(doubling_time) + '.png'

plt.savefig(figname, dpi=199)

plt.show()

print('Estimated doubling time: ' + str(doubling_time) + 'h')

print()

print('Linear fit results:')
print(out.fit_report(min_correl=0.5))


#get area under curve (integral viable cell density): over time and in total
fig = plt.figure(4)
auc_arr = np.zeros(hours.shape[0])
for i in range(hours.shape[0]): #i begins at 0
    auc_arr[i] = simps(np.array(viable[0:i+1]),hours[0:i+1])
plt.errorbar(hours, auc_arr, yerr=viable_err, fmt='-o')
plt.xlabel('days')
plt.ylabel('integral viable cell density\n(millions of cells*day/ml)')
figname = fignamebase + '_AUC.png'
plt.savefig(figname, dpi=199)
plt.show()
auc = simps(np.array(viable),hours)
print("Total integral viable cell density: " + '%.1f' % auc + " million cells*day/ml")


