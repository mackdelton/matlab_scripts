# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 14:50:07 2017

@author: lparker
"""

import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#Attempt to read in and visualize MAB data
data = pd.read_csv("./legacy_csv/testdata_20_100_2_rndgamma.csv", header=None,names=['id','selx','sely','rng','out'])
locs = pd.read_csv("./legacy_csv/testloc_20_100_2_rndgamma.csv", header=None, names=['alocx','alocy','blocx','blocy'])

N=len(locs) #Number of locations visited
tot = len(data) #Number of trials
data[['selx','sely']] = data[['selx','sely']].astype(float)
data[['id','out']] = data[['id','out']].astype(int)

#Indices are 1-based since being pulled from MATLAB

#Generate a new DataFrame containing the candidate locations and statistics of each
#candLoc = locs.loc[:,['alocx','alocy']]

#Concatenate candidate locations from 'locs' and their associated output(success and failures)
candLoc = locs.loc[:,['alocx','alocy']]
candLoc = pd.concat([(locs.loc[:,['alocx','alocy']]), pd.DataFrame({'out':[]}).astype(object), pd.DataFrame({'trials':[]}).astype(object), pd.DataFrame({'success':[]}).astype(object)], axis=1)

#candLoc.loc[[1],['out']] = np.array(data.loc[(data['id'] == 1),'out'])
tally = 0
for r in np.arange(0,N):
    candLoc.set_value(r, 'out', np.array(data.loc[data['id'] == (r+1),'out']))
    candLoc.set_value(r, 'trials', len(np.array(data.loc[data['id'] == (r+1),'out'])))
    candLoc.set_value(r, 'success', sum(np.array(data.loc[data['id'] == (r+1),'out'])))    
    #tally = tally + len(np.array(data.loc[data['id'] == (r+1),'out']))

#Histogram of frequency of visited locations
hdata = np.histogram(data['id'],bins=N)
hdata = hdata[0]
#bdata = np.arange(0,)

plt.close('all')

# Just a figure and one subplot
f, (ax1, ax2) = plt.subplots(2, 1)
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 10}
#Set up 2x1 plots for visualization
NN = 2
params = plt.gcf()
pltSize = params.get_size_inches()
params.set_size_inches( (pltSize[0]*NN, pltSize[1]*NN) )
#mpl.rc('font', **font)
mpl.rc('xtick', labelsize=10)
mpl.rc('ytick', labelsize=10) 

#Remove the extra space around the plot
plt.subplots_adjust(hspace=0.4)
#Plot LHS plot reflecting success and failures across locations
ind = np.arange(1,N+1)    # The x-axis for the locations considered
width = 0.35 #Width of bar plot
good = candLoc.loc[:,'success'].astype(int)
bad = candLoc.loc[:,'trials'].astype(int) - candLoc.loc[:,'success'].astype(int)
p1 = ax1.bar(ind, bad, width, color='r', align='center')
p2 = ax1.bar(ind, good, width, color='b', align='center')
ax1.set_xlabel('Candidate Receive Locations', multialignment='center', fontsize=14)
ax1.set_ylabel('Transmissions', multialignment='center', fontsize=14)
ax1.set_title('Projected Acomms Success/Failures', fontsize=14, fontweight='bold')
ax1.set_xticks(ind)
ax1.set_yticks(np.arange(0, 12, 2))
ax1.legend((p1[0], p2[0]), ('Failure', 'Success'))

#Plot RHS plot reflecting spatial distribution of trials across locations
colors = candLoc['success']/tot #Scale by 
#plt.scatter(locs['alocx'], locs['alocy'], s=tot*hdata, c=colors, alpha=0.5)
ax2.scatter(locs['alocx'], locs['alocy'], s=tot*hdata, alpha=0.5)
ax2.set_xlabel('Easting [m]', fontsize=14)
ax2.set_ylabel('Northing [m]', fontsize=14)
ax2.set_title('Spatial Distribution of Visits to Candidate Locations', fontsize=14, fontweight='bold')
ax2.set_xlim(min(locs['alocx'])-2, max(locs['alocx'])+1)
#plt.colorbar()
#plt.clim(0,1)

f.savefig('test.png', bbox_inches='tight')  
f.show()

print("%% Success: %6.2f \n"% (100*candLoc.loc[:,'success'].astype(int).sum(axis=0)/tot))
#plt.show()
