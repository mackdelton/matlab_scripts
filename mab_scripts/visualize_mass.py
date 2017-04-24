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


def parseMabInput(paraIn):
    #Extract relevant dataset pertaining to test scenario
    #Try to concatenate portions of dataIn relevant to the scenario of interest
    dataParse = dataIn.loc[ (dataIn['bernCnt'] == testIn[0]) & (dataIn['distribution'] == testIn[3]) &
                            (dataIn['soc'] == testIn[1]) & (dataIn['loc_cnt'] == testIn[4]) & 
                            (dataIn['iter'] == testIn[2]) & (dataIn['sol_type'] == testIn[5]),colT].reset_index(drop=True)
    #Generate new dataframe that stores statistics of each candidate location.
                            #Initialize dataframe with candidate locations and append 1) the resulting output vector (variable length)
                            #2) the total number of trials and 3) the number of successes.
    candLoc = pd.concat([(locs.loc[:,['alocx','alocy']]), pd.DataFrame({'out':[]}).astype(object), pd.DataFrame({'trials':[]}).astype(object), pd.DataFrame({'success':[]}).astype(object)], axis=1)
    
    for r in np.arange(0,N):
        #Extract the output for all instances of the "r"-th location
        candLoc.set_value(r, 'out', np.array(dataParse.loc[dataParse['id'] == (r+1),'out']))
        #Calculate the total length (or trials) associated with the "r"-th location
        candLoc.set_value(r, 'trials', len(np.array(dataParse.loc[dataParse['id'] == (r+1),'out'])))
        #Calculate the total number of successes at the "r"-th location    
        candLoc.set_value(r, 'success', sum(np.array(dataParse.loc[dataParse['id'] == (r+1),'out'])))    
    
        #Redundant command left here for posterity or in case it is useful later
        #tally = tally + len(np.array(data.loc[data['id'] == (r+1),'out']))

    return (candLoc, dataParse);

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
plt.close('all')

#STEP 0: Set up parameters of interest

#Set flag to determine whether you are collecting stats or just plotting examples

#STEP 0.1: If plotFeature == 1, re-define all parameters below as scalars
#(i.e., remove all tuples)
plotFeature = 0 #Set this value to 1
#Initialize parameters
#Define parameters of interest
#******************************************************************************
#bC = 1
bC = (1,2,3)       #Number of Bernoulli trials considered (1 - 1of1, 2 - 1of5, 3 - 1of10)
#*****************
soc = 0 #Type of Success of Communication curve (0 - Gamma, 1 - Exponential)
#soc = (0, 1)
#*****************
#nIter = 100  #Number of iterations (timesteps or duration, 100, 500, or 1000)
nIter = (100, 500, 1000)
#*****************
soln = 1     #Solution of interest (Gittins Index, GI - 0, Uninformed Random, UR - 1; Educated Guess, EG - 2)
#soln = (0, 1, 2)
#*****************
distrT = 0    #Distribution type (Uniform Random - 0, Mesh grid - 1)
#distrT = (0, 1)
#*****************
bStationary = 0 #Flag indicating whether or not the TX agent ("B") is stationary (1) or not ()
#bStationary = (0, 1)
#*****************
N = 20       #Number of candidate locations
#N = (20, 100)
#******************************************************************************
#Import test scenario file
# (e.g., 20 randomly distributed candidate locations: "*_rnd20.txt") and

if distrT == 0: #If Uniformly randomly distributed "arm" locations
    scen = 'rnd' #Define the relevant string
else:
    scen = 'mesh'

locsFileIn = './testloc_'+scen+str(N)+'_cmp.txt'
locs = pd.read_table(locsFileIn, header=None, sep=',', names=['alocx','alocy','blocx','blocy'])
#Update to match candidate location input file

#Create string for input file
dataFileIn = './improv_testdata10APR17_stationaryB_'+str(bStationary)+'.txt'
colT = ['id','selx','sely','rngP','rng','out','bernCnt','distribution','soc','loc_cnt','iter','sol_type','dist_tot'] #Define the column headers for data extraction
dataIn = pd.read_table(dataFileIn, header=None, sep=',', index_col=False, names=colT)
#Change data types of specific columns (i.e., float for some columns and int for others)
dataIn[['selx','sely','rngP','rng']] = dataIn[['selx','sely','rngP','rng']].astype(float)
dataIn[['id','out','bernCnt','distribution','soc','loc_cnt','iter','sol_type','dist_tot']] = dataIn[['id','out','bernCnt','distribution','soc','loc_cnt','iter','sol_type','dist_tot']].astype(int)

if plotFeature == 0:
    #Testing different sets of conditions
    #for b in bC:
    for b in bC:
        for nI in nIter:
            #Define tuple of inputs to parsing function
            testIn = (b, soc, nI, distrT, N, soln)
            locsStat, dataStat = parseMabInput(testIn)
            successTotal = locsStat.loc[:,'success'].astype(int).sum(axis=0)
            trialTotal = locsStat.loc[:,'trials'].astype(int).sum(axis=0)
            print("%% Success: %6.2f --- Distance: %7.2f\n"% ((100*successTotal/trialTotal), dataStat.loc[1,'dist_tot'] ))
else:
#****************** PLOTTING *********************************************

    #print(locsStat)
    #bdata = np.arange(0,)
    #plt.bar(hdata[0],hdata[1]) ????
    
    #Histogram of frequency of visited locations
    hdata = np.histogram(dataStat['id'],bins=N)
    hdata = hdata[0]
    
    # Just a figure and one subplot
    f, (ax1, ax2, ax3) = plt.subplots(3, 1)
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
    
    scl = 1   #Scale factor for plotting (1 - nIter=100, 0.05 - nIter=500, 0.001 - nIter=1000)
    #Remove the extra space around the plot
    plt.subplots_adjust(hspace=0.4)
    #Plot LHS plot reflecting success and failures across locations
    ind = np.arange(1,N+1)    # The x-axis for the locations considered
    width = 0.35 #Width of bar plot
    good = locsStat.loc[:,'success'].astype(int)
    bad = locsStat.loc[:,'trials'].astype(int) - locsStat.loc[:,'success'].astype(int)
    p2 = ax1.bar(ind, good, width, color='b', align='center')
    p3 = ax1.bar(ind, bad, width, bottom=good, color='r', align='center')
    ax1.set_xlabel('Candidate Receive Locations', multialignment='center', fontsize=14)
    ax1.set_ylabel('Transmissions', multialignment='center', fontsize=14)
    ax1.set_title('Projected Acomms Success/Failures', fontsize=14, fontweight='bold')
    ax1.set_xticks(ind)
    ax1.set_yticks(np.arange(0, max(locsStat['success']), max(locsStat['success'])/scl))
    #ax1.set_yticks(np.arange(0, 600, 10))
    ax1.legend((p2[0], p3[0]), ('Success','Failure'))
    
    #Plot RHS plot reflecting spatial distribution of trials across locations
    colors = locsStat['success']/nIter #Scale by 
    #plt.scatter(locs['alocx'], locs['alocy'], s=nIter*hdata, c=colors, alpha=0.5)
    #Scale size of scattered data
    ax2.scatter(locs['alocx'], locs['alocy'], s=nIter*hdata*scl, alpha=0.5)
    ax2.set_xlabel('Easting [m]', fontsize=14)
    ax2.set_ylabel('Northing [m]', fontsize=14)
    ax2.set_title('Spatial Distribution of Visits to Candidate Locations', fontsize=14, fontweight='bold')
    ax2.set_xlim(min(locs['alocx'])-2, max(locs['alocx'])+1)
    #plt.colorbar()
    #plt.clim(0,1)
    #Histogram of ranges calculated between agentA and agent B
    ax3 = plt.hist(dataStat['rng'],bins=60)
    
    f.savefig('CMP_27MAR17_PKR.png', bbox_inches='tight')  
    f.show()
    #plt.show()
