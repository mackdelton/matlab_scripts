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

'''
Data Available for analysis (as of 12FEB17)
LOCATION FILES
testloc_rnd20.txt
testloc_rnd100.txt
testloc_mesh20.txt
testloc_mesh100.txt

DATA FILES
testdata12FEB17.txt
'''

#STEP 1: Import test scenario (e.g., 20 randomly distributed candidate locations: "*_rnd20.txt") and specify corresponding parameters
locs = pd.read_table("./testloc_rnd20.txt", header=None, sep=',', names=['alocx','alocy','blocx','blocy'])
#Update to match candidate location input file
N = 20       #Number of candidate locations
distT = 0    #Distribution type (Uniform Random - 0, Mesh grid - 1)

#Define parameters of interest
bC = 2       #Number of Bernoulli trials considered (1 - 1of1, 2 - 1of5, 3 - 1of10)
soc_type = 1 #Type of Success of Communication curve (0 - Gamma, 1 - Exponential)
nIter = 100  #Number of iterations (timesteps or duration, 100, 500, or 1000)
scl = 1   #Scale factor for plotting (1 - 100, 0.05 - 500, 0.001 - 1000)
soln = 0     #Solution of interest (Gittins Index, GI - 0, Uninformed Random, UR - 1; Educated Guess, EG - 2)


#STEP 2: Read in and visualize MAB data
colT = ['id','selx','sely','rng','out','bernCnt','distribution','soc','loc_cnt','iter','sol_type'] #Define the column headers for data extraction
dataIn = pd.read_table("./testdata12FEB17.txt", header=None, sep=',', index_col=False, names=colT)
#Change data types of specific columns
dataIn[['selx','sely','rng']] = dataIn[['selx','sely','rng']].astype(float)
dataIn[['id','out','bernCnt','distribution','soc','loc_cnt','iter','sol_type']] = dataIn[['id','out','bernCnt','distribution','soc','loc_cnt','iter','sol_type']].astype(int)

#STEP 3: Extract relevant dataset pertaining to test scenario
#Try to concatenate portions of dataIn relevant to the scenario of interest
dataParse = dataIn.loc[ (dataIn['bernCnt'] == bC) & (dataIn['distribution'] == distT) &
                        (dataIn['soc'] == soc_type) & (dataIn['loc_cnt'] == N) & 
                        (dataIn['iter'] == nIter) & (dataIn['sol_type'] == soln),colT].reset_index(drop=True)

candLoc = pd.concat([(locs.loc[:,['alocx','alocy']]), pd.DataFrame({'out':[]}).astype(object), pd.DataFrame({'trials':[]}).astype(object), pd.DataFrame({'success':[]}).astype(object)], axis=1)

#candLoc.loc[[1],['out']] = np.array(data.loc[(data['id'] == 1),'out'])
tally = 0
for r in np.arange(0,N):
    #Extract the output for all instances of the "r"-th location
    candLoc.set_value(r, 'out', np.array(dataParse.loc[dataParse['id'] == (r+1),'out']))
    #Calculate the total length (or trials) associated with the "r"-th location
    candLoc.set_value(r, 'trials', len(np.array(dataParse.loc[dataParse['id'] == (r+1),'out'])))
    #Calculate the total number of successes at the "r"-th location    
    candLoc.set_value(r, 'success', sum(np.array(dataParse.loc[dataParse['id'] == (r+1),'out'])))    

    #Redundant command left here for posterity or in case it is useful later
    #tally = tally + len(np.array(data.loc[data['id'] == (r+1),'out']))

#Histogram of frequency of visited locations
hdata = np.histogram(dataParse['id'],bins=N)
hdata = hdata[0]
#bdata = np.arange(0,)

plt.close('all')
print("%% Success: %6.2f \n"% (100*candLoc.loc[:,'success'].astype(int).sum(axis=0)/nIter))
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
ax1.set_yticks(np.arange(0, max(candLoc['success']), max(candLoc['success'])/scl))
#ax1.set_yticks(np.arange(0, 600, 10))
ax1.legend((p1[0], p2[0]), ('Failure', 'Success'))

#Plot RHS plot reflecting spatial distribution of trials across locations
colors = candLoc['success']/nIter #Scale by 
#plt.scatter(locs['alocx'], locs['alocy'], s=nIter*hdata, c=colors, alpha=0.5)
#Scale size of scattered data
ax2.scatter(locs['alocx'], locs['alocy'], s=nIter*hdata*scl, alpha=0.5)
ax2.set_xlabel('Easting [m]', fontsize=14)
ax2.set_ylabel('Northing [m]', fontsize=14)
ax2.set_title('Spatial Distribution of Visits to Candidate Locations', fontsize=14, fontweight='bold')
ax2.set_xlim(min(locs['alocx'])-2, max(locs['alocx'])+1)
#plt.colorbar()
#plt.clim(0,1)

f.savefig('test.png', bbox_inches='tight')  
f.show()
#plt.show()
