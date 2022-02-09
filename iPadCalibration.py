# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:11:40 2021

@author: rithv
"""
#import sys
import daq     # interface to NI DAQ
import time
#import notify
import pickle
#import logging
import scipy.io
import scipy.signal
#import threading
import numpy             as     np
import scipy             as     sp
import datetime          as     dt
#import sounddevice       as     sd
import matplotlib.pyplot as     plt

def save():
	# save data to file
	print('Saving results to "%s"' % dataFile)
	params = ['repeats', 'outFs', 'inpFs', 'syncDur', 'syncAmp', 'sync', 'stim',
		'outSig', 'inpSum', 'inpSumFilt', 'inpSumSync', 'inpSumTrim',
		'inpAvg', 'startTime', 'stopTime', 'outMode', 'stimSL', 'stimSPL',
		'count', 'inpAll', 'cal0dB', 'thrSPL', 'gapDur', 'syncSPL', 
		'inpAvgTrim', 'inpAvgFilt', 'inpAvgSync']
	data = {}
	for param in params:
		if param in globals().keys():
			data[param] = globals()[param]

	pickle.dump(data, open(dataFile + '.P', 'wb'))
	sp.io.savemat(dataFile + '.mat', data)

time.sleep(5)

repeats = 100
count = 0

stimSPL = 50 #Must change for every SPL.
cal0dB = 103

if stimSPL > 55:
    syncSPL = 75
    inpRange = (-2,2)
elif stimSPL > 20:
    syncSPL = 65
    inpRange = (-.35,.35)
else:
    syncSPL = 55
    inpRange = (-.2,.2)
    
#Starting Procedure
#Time
print('Starting calibration at %.1f dB SPL' % stimSPL)
timeFormat = '%Y%m%d-%H%M%S'
startTime  = dt.datetime.now().strftime(timeFormat)
tic        = dt.datetime.now()
# Calibration
dataFile = 'filteredData/calibration-%s-%.1fdbspl%d' % (startTime, stimSPL, repeats)
inpFs    = 50e3       # input sampling frequency (Hz)
stimDur  = 5.243      # stim signal duration (s)
syncDur  = 10e-3      # sync signal duration (s)
gapDur   = 100e-3     # gap duration between sync and stim signals (s)
syncAmp  = 10**((syncSPL-cal0dB)/20) if syncSPL else 0  # sync amplitude (V)


# Creating the sync for reference
t        = np.arange(0, syncDur, 1/inpFs)
sync     = np.exp(-(t-syncDur*2/4)**2/(2*(syncDur/50)**2)) * syncAmp
sync    -= np.exp(-(t-syncDur*3/4)**2/(2*(syncDur/50)**2)) * syncAmp
sync    -= np.exp(-(t-syncDur*1/4)**2/(2*(syncDur/50)**2)) * syncAmp
    
inpTotalN  = int(((syncDur + gapDur + stimDur + gapDur) * repeats) * inpFs + 1)
inpMLSN    = int(stimDur * inpFs)
inpSum     = np.zeros((2, inpTotalN))     # original recordings
inpSumFilt = np.zeros((2, inpTotalN))     # bandpass filtered
inpSumSync = np.zeros((2, inpMLSN))     # including start of sync frame
inpSumTrim = np.zeros((2, inpMLSN))     # excluding sync frame and gap


if 'inpAll' in locals(): del inpAll
    
################################################################################
print("Starting to read DAQ input")

#Left Ear: Channel 2
#Right Ear: Channel 1

ai = daq.AnalogInput ('/dev1/ai0:1', inpFs, inpTotalN, range=inpRange, accurateFs=True)

try:
    # read the recorded signal
    ai.start()
    #Reads the entire recording that is being played. 
    inpSig = ai.read()
    
    print("Starting Signal Processing")
    #Start Signal processing
    inpSum[:, :inpSig.shape[-1]] += inpSig
    
    # filter
    b, a = sp.signal.butter(3, np.array([20,20e3])/(inpFs/2), 'bandpass')
    inpFiltMLS = sp.signal.filtfilt(b, a, inpSig) # for extracting MLS
    
    b, a = sp.signal.butter(3, np.array([1,1000])/(inpFs/2), 'bandpass')
    inpFiltSync = sp.signal.filtfilt(b, a, inpSig) # for correlating with sync frame
    #inpSumFilt[:, :inpFilt.shape[-1]] += inpFilt
    
    #CrossCorrelation
    
    
    xcor = sp.signal.correlate(inpFiltSync[0,:], sync, 'same') 
    maxCorr = np.max(xcor)
    lower = .95*maxCorr
    upper = 1.05*maxCorr
    d = int((syncDur + gapDur + stimDur + gapDur) * inpFs)
    peaks,_ = sp.signal.find_peaks(xcor,height = (lower,upper), distance = d)
     
    for i in peaks:
        #Indecies for the Sync Frame
        indStart = int(i - int(inpFs*syncDur/2))
        indEnd = int(i +  int(inpFs*syncDur/2))
        #Indecies for the MLS
        indGS = int(indEnd + int(inpFs*gapDur)) # starts at the end of gap
        indGE = int(indGS + int(inpFs*stimDur))# ends at the gap after the MLS
        
        inpSumSync = inpFiltMLS[:,indGS + 1:indGE + 1] #Extracts the data between the two gaps
        inpSumTrim[:,:inpSumSync.shape[-1]] += inpSumSync
        count = count + 1;
       
    # calculate average of accumulated input signals
    inpAvg     = inpSum     / count
    inpAvgFilt = inpSumFilt / count
    inpAvgSync = inpSumSync / count
    inpAvgTrim = inpSumTrim / count 
    stopTime = dt.datetime.now().strftime(timeFormat)
    save()
    toc = dt.datetime.now()
    print('Total time: %.2f s' % (toc-tic).total_seconds())
    
    ################################################################################
    print("Starting Plotting")
    # plotting
    fig, (ax1,ax2,ax3) = plt.subplots(3,1, sharex=True)#, sharey=True)
    fig.set_size_inches(12, 8)
    
    #Original Signal Plot
    sig = inpFiltMLS
    t = np.arange(0, sig.shape[-1]) / inpFs
    ax1.plot(t, sig[0,:])
    ax1.set_title('Input Signal with Sync Frame Left Ear')
    ax1.set_ylabel('Voltage (V)')
    ax1.grid(True)
    
    ax2.plot(t, sig[1,:])   
    ax2.set_title('Input Signal with Sync Frame Right Ear')
    ax2.set_ylabel('Voltage (V)')
    ax2.grid(True)
    
    #Extracted Signal Average Plot
    sig = inpAvgTrim
    t = np.arange(0, sig.shape[-1]) / inpFs
    for i in range(sig.shape[0]):
    		ax3.plot(t, sig[i,:])
    ax3.set_title('Average of MLS')
    ax3.set_ylabel('Voltage (V)')
    ax3.grid(True)
    
    plt.show()
    
except:
    raise
    
finally:
    del ai
    
    

