# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:12:51 2020

@author: miseminger
"""

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



'''these lines should be filled in before running the code'''

filename = r'D:\\F8Acti2l.xlsx' #name of Excel file to read from
sheet = 0 #sheet name if applicable (0 is default)
start_row = 2 - 2 #number of the starting row in your Excel file, minus 2
end_row = 10 - 2 #number of the final row in your Excel file, minus 2
linearpoints = 3 #number of points in linear phase to fit
first_linearpoint = 0 #first point is 0
title = 'F8Acti2l' #plot title
errorbars = False #set to True if you have multiple samples per data point 
plot_NOVA = True #set to True if you want plots of glucose, lactate, BUN and pH; otherwise set to False
feedmarkers = True
samplemarkers = True





'''now the code begins'''

#set basis for figure names
#fignamebase = '/Users/miseminger/Desktop/graphs/' + title + '_' + str(start_row) + '_' + str(end_row)

fignamebase = r'D:\\' + title + '_'

#read in the data from Excel!
skip_rows = start_row - 1 #number of rows to skip
df = pd.read_excel(filename, sheet_name=sheet)
days = df['Date finished'][start_row:end_row+1]
viability = df['Viability'][start_row:end_row+1]  #percent viability
viable = df['Viable Cell Conc.'][start_row:end_row+1]/1000000 #millions of cells per ml
    

if errorbars == False:
    viability_err = 0
    viable_err = 0

else:
    viability_err = df['Mean_Viability'][start_row:end_row+1]
    viable_err = df['Mean_Viable'][start_row:end_row+1]/1000000

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


#plot % viability

fig = plt.figure(1)

plt.errorbar(hours,viability, yerr=viability_err, color='mediumspringgreen', fmt='-o')       
plt.xlabel('days')
plt.ylabel('% viable')
plt.ylim(0,110)


#show sampling and feeding days with triangles on viability plot

if feedmarkers==True:
    feed = df['feed'][start_row:end_row+1] * (viability + 5) #feeding days
    plt.scatter(hours, feed, color='mediumaquamarine', marker='v', label='added feed')

if samplemarkers==True:
        sample = df['sample'][start_row:end_row+1] * (viability - 5) #sampling days
        plt.scatter(hours, sample, color='aquamarine', marker='^', label='sampled CM')

plt.legend(loc='lower left')

figname = fignamebase + '_percentviable.png'
plt.savefig(figname, dpi=199)


#plot viable cell concentration in millions/ml

fig = plt.figure(2)

plt.errorbar(hours,viable, yerr=viable_err, fmt='-o')
plt.xlabel('days')
plt.ylabel('viable cells (millions/ml)')

#show sampling days with black triangles on viable cell concentration plot

if feedmarkers==True:
    feed = df['feed'][start_row:end_row+1] * (viable + 0.3) #feeding days
    plt.scatter(hours, feed, color='mediumaquamarine', marker='v', label='added feed')

if samplemarkers==True:
        sample = df['sample'][start_row:end_row+1] * (viable - 0.3) #sampling days
        plt.scatter(hours, sample, color='aquamarine', marker='^', label='sampled CM')

plt.legend(loc='lower right')

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
y = slope * np.arange(-1,10) + intercept  
doubling_time = int(np.log(2)/slope*24) #doubling time in hours

plt.plot(np.arange(-1,10), y, 'r-', label='linear fit')
plt.errorbar(hours, viable_log, yerr=viable_err/viable, label='viable cells', fmt='-o')  
plt.xlabel('days')
plt.ylabel('ln(millions of cells/ml)')
plt.legend(loc='upper left')
plt.ylim(np.min(viable_log)-1, np.max(viable_log)+1)
plt.xlim(-0.7,)
plt.text(7, np.min(viable_log), 'estimated doubling time: ' + str(doubling_time) + 'h', {'color': 'k', 'fontsize': 10})

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

plt.errorbar(hours, auc_arr, yerr=viable_err, fmt='-o', color='indigo')
plt.xlabel('days')
plt.ylabel('integral viable cell density\n(millions of cells*day/ml)')

figname = fignamebase + '_AUC.png'
plt.savefig(figname, dpi=199)

plt.show()

auc = simps(np.array(viable),hours)
print("Total integral viable cell density: " + '%.1f' % auc + " million cells*day/ml")





#plot NOVA data if applicable



if plot_NOVA == True:

    glucose = df['glucose'][start_row:end_row+1] * 0.18016 #glucose in g/l
    lactate = df['lactate'][start_row:end_row+1] * 0.08907 #lactate in g/L
    BUN = df['BUN'][start_row:end_row+1]
    pH = df['pH'][start_row:end_row+1]


    #plot glucose and lactate in g/l

    fig = plt.figure(5)

    plt.scatter(hours, glucose, label='glucose', color='darkcyan')
    plt.scatter(hours, lactate, label='lactate', color='lightcoral')
    plt.ylabel('concentration (g/L)')
    plt.xlabel('days')
    plt.legend()
    plt.grid()

    figname = fignamebase + '_glucose_lactate.png'
    plt.savefig(figname, dpi=199)

    plt.show()

    

    #plot ammonium in mM

    fig = plt.figure(6)

    plt.scatter(hours, BUN, color='orange')
    plt.ylabel('BUN concentration (mM)')
    plt.xlabel('days')
    
    figname = fignamebase + '_BUN.png'
    plt.savefig(figname, dpi=199)

    plt.show()



    #plot pH

    fig = plt.figure(7)

    plt.scatter(hours, pH, color='darkturquoise')
    plt.ylabel('pH')
    plt.xlabel('days')
    plt.grid()

    figname = fignamebase + '_pH.png'
    plt.savefig(figname, dpi=199)

    plt.show()
    
    
    #plot target protein concentration vs. integral cell concentration
    
    fig = plt.figure(8)
    
    prot_nM = [4.513230047, 19.80373578, 55.46120443, 99.65672088, 133.5302892, 145.1529999, 161.8156372, 93.17130731]

    plt.errorbar(auc_arr[1:], prot_nM, yerr=0, fmt='-o', color='r')   
    plt.ylabel('[prot] (nM)')
    plt.xlabel('integral viable cell concentration \n(millions of cells/ml*day)')
    plt.grid()

    figname = fignamebase + '_proteinvscells.png'
    plt.savefig(figname, dpi=199)

    plt.show()
    
    
        #plot protein vs. time
    
    fig = plt.figure(9)
    
    prot_nM = [4.513230047, 19.80373578, 55.46120443, 99.65672088, 133.5302892, 145.1529999, 161.8156372, 93.17130731]

    plt.errorbar(hours[1:], prot_nM, yerr=0, fmt='-o', color='k')   
    plt.ylabel('[prot] (nM)')
    plt.xlabel('days')
    plt.grid()

    figname = fignamebase + '_proteintimecourse.png'
    plt.savefig(figname, dpi=199)

    plt.show()
    
    
        #plot protein in g/L vs. integral cell concentration
    
    fig = plt.figure(10)
    
    prot_g = np.array(prot_nM) * 0.168
    
    plt.errorbar(auc_arr[1:], prot_g, yerr=0, fmt='-o', color='r')   
    plt.ylabel('[prot] (mg/L)')
    plt.xlabel('integral viable cell concentration \n(millions of cells/ml*day)')
    plt.grid()

    figname = fignamebase + '_proteinvscells_mg.png'
    plt.savefig(figname, dpi=199)

    plt.show()
    
    
        #plot protein in g/L vs. time
    
    fig = plt.figure(11)
    
    plt.errorbar(hours[1:], prot_g, yerr=0, fmt='-o', color='k')   
    plt.ylabel('[prot] (mg/L)')
    plt.xlabel('days')
    plt.grid()

    figname = fignamebase + '_proteintimecourse_mg.png'
    plt.savefig(figname, dpi=199)

    plt.show()
