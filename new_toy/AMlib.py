import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from scipy import special
import scipy.optimize as sco
import scipy.signal as scs
from glob import glob 
import os 
import re
import pickle as pkl

def expon1(t,tau_rise=20,tau_fall=50,amp=1,t_offset=0):
    #return expon3(t,tau_rise,tau_fall,1,1,amp,0,0,t_offset)
    # Don't call expon3 to if not using extra fallt times, speeds up execution by 50%
    exp = np.zeros(t.shape[-1])
    if t_offset>t[-1]:
        return exp

    t2 = t[t>t_offset]
    exp[t>t_offset] = amp*np.exp(-(t2-t_offset)/(tau_fall)) - (amp)*np.exp(-(t2-t_offset)/tau_rise)

    return exp

def fitPulse(t,trace,p0=None,tau_fall=None,tau_rise=None,printDebug=False,returnFunc=False):
    '''
    Algorithm to fit an exponential pulse to a trace. Can provide a fixed tau_rise or tau_fall, both or neither.
    If p0 is not provided, it will make initial guess based on where trace crosses 50% of max
    printDebug enabled will print this guess

    Returns: fit, cov, fitNames, (if returnFunc, the expon function used)
    '''

    tfQ = tau_fall is None # Should fall be fit?
    trQ = tau_rise is None # Should rise be fit?

    if trace.ndim!=1:
        raise ValueError('Only 1D trace supported')
    if len(t)!=len(trace):
        raise ValueError('t and trace must be same shape')

    if p0 is None:
        #Guess taus based on crossing 50% of max
        amax = np.argmax(trace)
        maxV = np.max(trace)
        r50 = np.arange(amax)[trace[:amax][::-1]<maxV/2][0]  # time from 50% max to max
        f50 = np.arange(len(trace)-amax)[trace[amax:]<maxV/2][0]  #time from max to 50% max
        guessRise = t[f50]
        guessFall = t[r50]

        guessAmp=maxV*3  #There's an actual scaling for this, might be good to implement
        guessStart = amax-2*r50 #Left as index, converted to time in p0
        guessBase = np.median(trace[:guessStart-5*r50])

    #Modify the fitting functions to use delta_tau_fall such that tau_fall>tau_rise, but only when fitting tau_fall.
    if (tfQ&trQ):  #Neither provided, fit both
        def exponI(t,tau_rise,dtau_fall,amp,t_start,baseline):
            tau_fall=tau_rise+dtau_fall
            return expon1(t,tau_rise,tau_fall,amp,t_start)+baseline
        p0=[guessRise,max([guessFall-guessRise,guessRise/10]),guessAmp,t[guessStart],guessBase]
        fitNames = ['tau_rise','tau_fall','amp','t_offset','baseline']
        def exponI2(t,tau_rise,tau_fall,amp,t_start,baseline):  #Return function
            return expon1(t,tau_rise,tau_fall,amp,t_start)+baseline
    elif trQ:  #Fixed fall, fit rise
        def exponI(t,tau_rise,amp,t_start,baseline):
            return expon1(t,tau_rise,tau_fall,amp,t_start)+baseline
        p0=[guessRise,guessAmp,t[guessStart],guessBase]
        fitNames = ['tau_rise','amp','t_offset','baseline']
        exponI2 = exponI
    elif tfQ:  #Fixed rise, fit fall
        def exponI(t,dtau_fall,amp,t_start,baseline):
            tau_fall=tau_rise+dtau_fall
            return expon1(t,tau_rise,tau_fall,amp,t_start)+baseline
        p0=[max([guessFall-tau_rise,tau_rise/10]),guessAmp,t[guessStart],guessBase]
        fitNames = ['tau_fall','amp','t_offset','baseline']
        def exponI2(t,tau_fall,amp,t_start,baseline):
            return expon1(t,tau_rise,tau_fall,amp,t_start)+baseline
    else: #Both fixed
        def exponI(t,amp,t_start,baseline):
            return expon1(t,tau_rise,tau_fall,amp,t_start)+baseline
        p0=[guessAmp,t[guessStart],guessBase]
        fitNames = ['amp','t_offset','baseline']
        exponI2=exponI

    if printDebug:
        print(p0)
    bounds = np.array([[0,np.inf]]*(len(p0)-1) + [[-np.inf,np.inf]]).T

    try:
        fit,cov = sco.curve_fit(exponI,t,trace,p0=p0,bounds=bounds,maxfev=20000)
    except:
        print(p0)
        raise
    if 'tau_fall' in fitNames:
        if 'tau_rise' in fitNames:
            tau_rise = fit[fitNames.index('tau_rise')]
        fit[fitNames.index('tau_fall')]+=tau_rise  #Change from dtau_fall to tau_fall

    if returnFunc:
        return fit,cov,fitNames,exponI2
    else:
        return fit,cov,fitNames


def alignPulses(pulse,alignTo='self',p=0.2,thresh=50,maxShift=100,returnOnlyRoll=False,alignType='pMax',printRoll=False):
    #alignType can be  pMax, max, threshold, input
    #if input is selected, a value/array must be input under alignTo which matches the pulse dimensions.

    if alignType not in ['pMax', 'max', 'threshold','input']:
        print('Unrecognized alignType. Please choose from pMax, max, threshold, input.' )
        return pulse

    if alignType == 'input':
        if pulse.ndim==1:
            if type(alignTo)==int:
                rollBy = alignTo
            if len(alignTo)==1:
                rollBy = alignTo[0]
            else:
                print('For single pulse, alignTo must be single integer value.')
                return pulse

            return np.roll(pulse,rollBy)
        else:
            if len(alignTo)!=len(pulse):
                print('alignTo must be same length has pulse array.')
                return pulse

            return np.array([np.roll(pulseI,rollBy) for pulseI,rollBy in zip(pulse,alignTo)])

    if alignTo=='self':
        if pulse.ndim==1:
            print('Align to self requires >1 dimension')
            return pulse
        align1 = alignPulses(pulse,alignTo=0,p=p,thresh=thresh,maxShift=pulse.shape[-1],returnOnlyRoll=True,alignType=alignType)
        alignTo = -int(np.median(align1))

    if pulse.ndim>1:
        return np.array([alignPulses(pulseI,alignTo,p,thresh,maxShift,returnOnlyRoll,alignType,printRoll=False) for pulseI in pulse])

    if alignType=='pMax':
        minIndex = np.max([0,alignTo-maxShift])
        maxIndex = np.min([len(pulse),alignTo+maxShift])
        maxV = np.max(pulse[minIndex:maxIndex])
        aboveTh = np.arange(minIndex,maxIndex)[pulse[minIndex:maxIndex]>maxV*p]
        if len(aboveTh)>0:
            rollBy = alignTo - aboveTh[0]
        else:
            rollBy = 0

    elif alignType=='max':
        minIndex = np.max([0,alignTo-maxShift])
        maxIndex = np.min([len(pulse),alignTo+maxShift])
        rollBy = alignTo - (np.argmax(pulse[minIndex:maxIndex])+minIndex)

    elif alignType=='threshold':
        minIndex = np.max([0,alignTo-maxShift])
        maxIndex = np.min([len(pulse),alignTo+maxShift])
        aboveTh = np.arange(minIndex,maxIndex)[pulse[minIndex:maxIndex]>thresh]
        if len(aboveTh)>0:
            rollBy = alignTo - aboveTh[0]
        else:
            rollBy = 0

    if printRoll:
        print('Rolling By:',rollBy)
    if returnOnlyRoll:
        return rollBy
    else:
        return np.roll(pulse,rollBy)

    
def subtractBaselines(events,nbins=500):
    baselines = np.median(events[...,:nbins],axis=-1)
    return events - baselines.reshape(baselines.shape+(1,))

def makeTemplate(trace,locs,tlen=4096,fs=None,tau_rise=None,tau_fall=None,printDebug=True):
    '''
    NOTE: Trace must have positive going pulses to align to. Baseline subtraction is applied to each separated trace.
    Returns dictionary with Trace-Mean, TemplateM (as fit to Trace-mean), Template (Normalized and centered), fit and cov
    '''
    traces = np.asarray([trace[lI-tlen//2:lI+tlen//2] for lI in locs])
    traces = subtractBaselines(traces,nbins=tlen//4)
    offsets = alignPulses(traces,p=0.5,alignTo=tlen//2,returnOnlyRoll=True)  #Aligning to %max on the rise is more sharply defined than to the max of the trace
    traces = np.asarray([trace[lI-tlen//2:lI+tlen//2] for lI in locs-offsets])  #Rather than roll initial trace, just grab the trace again with the correct location
    traces = subtractBaselines(traces,nbins=tlen//4)

    traceM = np.median(traces,axis=0)
    traceM*=1/np.max(traceM)

    taxis = np.arange(tlen)/fs
    if False:
        #Find intial guess
        tmax,tmaxA = np.max(traceM),np.argmax(traceM)
        p50U = taxis[np.argmax(traceM)] - taxis[np.arange(tlen)[traceM>np.max(traceM)*0.5][0]]
        p50D = taxis[np.arange(tlen)[traceM>np.max(traceM)*0.5][-1]] - taxis[np.argmax(traceM)]

        guess = (p50U,p50D,1,taxis[tmaxA]-2*p50U)
        fit,cov = sco.curve_fit(expon1,taxis,traceM,p0=guess,maxfev=10000)
        #print(guess)
    else:
        fit,cov,fitNames,efunc = fitPulse(taxis,traceM,p0=None,tau_fall=tau_fall,tau_rise=tau_rise,printDebug=printDebug,returnFunc=True)
        #Fit is (rise_time),(fall_time),amp,start_time,baseline  ; with rise/fall times if not provided as fixed

    #This version has different cut process, first convolves the trace with the template, which can be a pulse template, or any filter kernel.
    
    tlen = len(template)
    
    if True: # Verify sign
        if np.max(trace)-np.median(trace) < -np.min(trace)+np.median(trace):
            #print('flipped')
            trace = trace*-1  #Assume that larger deviations negative than positive means trace is inverted
    if True: #Ensure max at middle to avoid time shift
        template = np.roll(template,tlen//2-np.argmax(template))
    if True: #AC Couple template, this affects amplitude scaling but this doesn't matter for this process.
        template = template - np.mean(template)
    traceC = np.correlate(trace,template,mode='same')
    
    traceC = traceC[tlen//2:-tlen//2]  #Cut poorly defined region
    traceCs = traceC[:len(traceC)//tlen*tlen].reshape(len(traceC)//tlen,tlen)  #initial reshape into tlen chunks
    
    stdVs = np.std(traceCs,axis=-1)
    stdVsS = np.sort(stdVs)
    cutSTD = stdVs < stdVsS[len(stdVsS)//10]*2  #Take 2x the 10% smallest std as initial cut
    print(cutSTD)
    if plotHist:
        ax = plt.subplots(1,2,figsize=(11,5))[1]
        ax[0].set(title='STD Histogram')
        ax[1].set(title='FFTSum Histogram')
        
        ax[0].hist(stdVs,bins=50) #,range=(0,1e-6))
        ax[0].hist(stdVs[cutSTD],bins=50) #,range=(0,1e-6))
    
    for ii in range(3):
        cutSTD = stdVs<2*np.median(stdVs[cutSTD])
        if plotHist:
            plt.hist(stdVs[cutSTD],bins=50) #,range=(0,1e-6))
            print(sum(cutSTD),np.mean(stdVs[cutSTD]))
    
    resB = np.median(stdVs[cutSTD])

    for ii in range(5):
        #Multiple passes to establish baseline resolution
        for multI in range(2,10):
            #Raise threshold until at least 5 noise regions are found
            locs,peaks = scs.find_peaks(traceC,height=2*resB,prominence=multI*resB)
            noiseT = getNoiseBetweenPeaks(traceC,locs,removeOutliers=False)
            if len(noiseT)>minNoiseTraces:
                break
        if printDebug&(multI>2):
            print("Warn: raised threshold to",multI)
        if printDebug:
            print(f'ResB: {resB} ; Number of peaks found: {len(locs)} ; Number of Noise Traces: {len(noiseT)}')
        if len(noiseT)==0:
            print(resB)
            plt.show()
            plt.plot(traceC)
            plt.axhline(3*resB,c='r')
            plt.show()
            raise RuntimeError(f"No noise traces recovered ; Number of 'peaks' = {len(locs)}")
            
        resB = np.median(np.std(noiseT,axis=-1))
    
    locs += tlen//2  #Correct for offset due to convolution
    noiseT = getNoiseBetweenPeaks(trace,locs,removeOutliers=False)
    if usePeriodogram:
        freqs,noises = scs.periodogram(noiseT,fs=fs)
        noises = noises**0.5
    else: 
        freqs = np.fft.rfftfreq(tlen,d=1/fs)
        noises = np.abs(np.fft.rfft(noiseT))


    
    noiseSum = np.sum(noises[:,1:],axis=-1)
    cutNS = np.full(len(noiseSum),True)
    if plotHist:
        ax[1].hist(noiseSum,bins=20,range=(0,1e-15))
    for ii in range(3):
        cutNS = noiseSum < np.median(noiseSum[cutNS])+2*np.std(noiseSum[cutNS])
        if plotHist:
            ax[1].hist(noiseSum[cutNS],bins=20,range=(0,1e-15))

    #print(len(cutNS))
    noiseFFT = np.mean(noises[cutNS]**2,axis=0)**0.5
    noiseFFT[0]=np.inf

    if plotHist:
        plt.show()

    returnX = {'Freqs':freqs}
    if usePeriodogram:  #Are we in PSD units or FFT units, name appropriately.
        returnX['NoisePSD']=noiseFFT
    else:
        returnX['NoiseFFT']=noiseFFT

    if returnSingles:
        returnX['Singles'] = noises[cutNS]
    if returnTraces:
        returnX['Traces'] = noiseT[cutNS]
    if returnPeaks:
        returnX['Locs'],returnX['Heights'] = locs,peaks['peak_heights']

    if makeTemplateQ:
        pcut = (peaks['peak_heights']>resB*templateThreshold[0]) & (peaks['peak_heights']<resB*templateThreshold[1])
        print('pcut:',pcut)
        if sum(pcut)==0:
            print('Unable to find pulses to generate template')
            returnX['TemplateD'] = None
        else:
            #            try:
            #if True:
            #plt.plot(trace[pcut])
            #plt.title('trace')
            #plt.show()
            templateD = makeTemplate(trace,locs[pcut],tlen,fs,tau_fall=tau_fall,tau_rise=tau_rise,printDebug=printDebug)
            templateD['NTraces'] = sum(pcut)
            returnX['TemplateD'] = templateD
            #except:
            ##else:
            #    print('Unable to fit')
            #    returnX['TemplateD'] = None
    
    return returnX
    
# Find all regions in data where 4096-bin traces have no peaks
def getNoiseBetweenPeaks(trace,peaks,nLeft=100,nRight=400,n=4096,removeOutliers=True):
    '''
    trace: continuous trace
    peaks: locations where peaks have been found
    nLeft: Number of bins left of each peak that should be discarded
    nRight: Same, right of each peak
    n=Length of each trace
    removeOutliers: Should a simple cut based on trace maxes be applied.

    returns traces
    '''
    noise=[]
    for p1,p2 in zip(peaks[:-1],peaks[1:]):
        nn = (p2-p1-(nLeft+nRight))//n  # Count number of noise traces between peaks, including buffer
        for ni in range(nn):
            tstart = p1+nRight+ni*n
            tend = p1+nRight+(ni+1)*n
            noise.append(trace[tstart:tend])
    noise = np.asarray(noise)
    if removeOutliers:
        for passI in range(2):
            maxes = np.max(noise,axis=-1)-np.median(noise,axis=-1)
            means = np.mean(noise,axis=-1)
            cuts = (maxes<np.mean(maxes)+3*np.std(maxes))&(means<np.mean(means)+3*np.std(means))
            noise=noise[cuts]

    return noise


def makeNoiseT(trace,template,fs,minNoiseTraces=30,usePeriodogram=False,plotHist=True,returnSingles=False,returnTraces=False,
               returnPeaks=False,templateThreshold=(7,100),makeTemplateQ=True,tau_rise=None,tau_fall=None,
               printDebug=False):
    #This version has different cut process, first convolves the trace with the template, which can be a pulse template, or any filter kernel.

    tlen = len(template)

    if True: # Verify sign
        if np.max(trace)-np.median(trace) < -np.min(trace)+np.median(trace):
            #print('flipped')
            trace = trace*-1  #Assume that larger deviations negative than positive means trace is inverted
    if True: #Ensure max at middle to avoid time shift
        template = np.roll(template,tlen//2-np.argmax(template))
    if True: #AC Couple template, this affects amplitude scaling but this doesn't matter for this process.
        template = template - np.mean(template)
    traceC = np.correlate(trace,template,mode='same')

    traceC = traceC[tlen//2:-tlen//2]  #Cut poorly defined region
    traceCs = traceC[:len(traceC)//tlen*tlen].reshape(len(traceC)//tlen,tlen)  #initial reshape into tlen chunks

    stdVs = np.std(traceCs,axis=-1)
    stdVsS = np.sort(stdVs)
    cutSTD = stdVs < stdVsS[len(stdVsS)//10]*2  #Take 2x the 10% smallest std as initial cut
    if plotHist:
        ax = plt.subplots(1,2,figsize=(11,5))[1]
        ax[0].set(title='STD Histogram')
        ax[1].set(title='FFTSum Histogram')

        ax[0].hist(stdVs, bins=50) #,range=(0,1e-6))
        ax[0].hist(stdVs[cutSTD], bins=50) #,range=(0,1e-6))

    for ii in range(3):
        cutSTD = stdVs<2*np.median(stdVs[cutSTD])
        if plotHist:
            plt.hist(stdVs[cutSTD] ,bins=50) #,range=(0,1e-6))
            print(sum(cutSTD),np.mean(stdVs[cutSTD]))
    if plotHist: plt.show()
            
    resB = np.median(stdVs[cutSTD])

    for ii in range(5):
        #Multiple passes to establish baseline resolution
        for multI in range(2,10):
            #Raise threshold until at least 5 noise regions are found
            locs,peaks = scs.find_peaks(traceC,height=2*resB,prominence=multI*resB)
            noiseT = getNoiseBetweenPeaks(traceC,locs,removeOutliers=False)
            if len(noiseT)>minNoiseTraces:
                break
        if printDebug&(multI>2):
            print("Warn: raised threshold to",multI)
        if printDebug:
            print(f'ResB: {resB} ; Number of peaks found: {len(locs)} ; Number of Noise Traces: {len(noiseT)}')
        if len(noiseT)==0:
            print(resB)
            plt.show()
            plt.plot(traceC)
            plt.axhline(3*resB,c='r')
            plt.show()
            raise RuntimeError(f"No noise traces recovered ; Number of 'peaks' = {len(locs)}")

        resB = np.median(np.std(noiseT,axis=-1))

    locs += tlen//2  #Correct for offset due to convolution
    noiseT = getNoiseBetweenPeaks(trace,locs,removeOutliers=False)

    #plt.clf()
    #plt.plot(noiseT)
    #plt.show()

    if usePeriodogram: # PSD: power spectral density
        freqs,noises = scs.periodogram(noiseT,fs=fs)
        noises = noises**0.5 # sqrt(PSD)
    else:
        freqs = np.fft.rfftfreq(tlen,d=1/fs) # frequency bins
        noises = np.abs(np.fft.rfft(noiseT)) # FFT of a real input



    noiseSum = np.sum(noises[:,1:],axis=-1)
    cutNS = np.full(len(noiseSum),True)
    if plotHist:
        ax[1].hist(noiseSum,bins=20,range=(0,1e-15))
    for ii in range(3):
        cutNS = noiseSum < np.median(noiseSum[cutNS])+2*np.std(noiseSum[cutNS])
        if plotHist:
            ax[1].hist(noiseSum[cutNS],bins=20,range=(0,1e-15))

    #print(len(cutNS))
    noiseFFT = np.mean(noises[cutNS]**2,axis=0)**0.5
    noiseFFT[0]=np.inf

    if plotHist:
        plt.show()

    returnX = {'Freqs':freqs}
    if usePeriodogram:  #Are we in PSD units or FFT units, name appropriately.
        returnX['NoisePSD']=noiseFFT
    else:
        returnX['NoiseFFT']=noiseFFT

    if returnSingles:
        returnX['Singles'] = noises[cutNS]
    if returnTraces:
        returnX['Traces'] = noiseT[cutNS]
    if returnPeaks:
        returnX['Locs'],returnX['Heights'] = locs,peaks['peak_heights']

    if makeTemplateQ:
        pcut = (peaks['peak_heights']>resB*templateThreshold[0]) & (peaks['peak_heights']<resB*templateThreshold[1])
        if sum(pcut)==0:
            print('Unable to find pulses to generate template')
            returnX['TemplateD'] = None
        else:
            try:
            #if True:
                templateD = makeTemplate(trace,locs[pcut],tau_fall=tau_fall,tau_rise=tau_rise,printDebug=printDebug)
                templateD['NTraces'] = sum(pcut)
                returnX['TemplateD'] = templateD
            except:
            #else:
                print('Unable to fit')
                returnX['TemplateD'] = None

    return returnX

