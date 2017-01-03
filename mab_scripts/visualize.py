# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 14:50:07 2017

@author: lparker
"""

import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#Attempt to read in and visualize MAB data
data = pd.read_csv("./testdata_20_100_1_meshexp.csv", header=None,names=['id','selx','sely','rng','out'])
locs = pd.read_csv("./testloc_20_100_1_meshexp.csv", header=None, names=['alocx','alocy','blocx','blocy'])

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

colors = candLoc['success']/tot #Scale by 

plt.scatter(locs['alocx'], locs['alocy'], s=tot*hdata, c=colors, alpha=0.5)
plt.colorbar()
plt.clim(0,1)
#plt.show()
