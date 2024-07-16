'''
 - Full toy simulation of a train of pulses.
 - Analysis of the train of pulses:
    * noise
    * peak finding
    * resolution

P. Franchini 7/2024
'''

import os
import sys
import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.signal as scs

# import Tsepelin code
exec(open("../mod_helium3.py").read())

# import QUEST-DMC code
exec(open("../mod_quest.py").read())

# import Adam M. library
sys.path.append('.')
import AMlib as afp

#===========================================================

# Input:
time=1677021898
wire=1 # (1,2)

energy_pdf = np.arange(0, 5.0e6, 1)  # energy spectrum [ev]

rate = 0.01      # [events/second], rate of events
max_time = 3600  # [second], total lenght of the sample
sampling = 100    # [Hz], sampling (points per second)

filename = "output.txt"

verbose = False

#===========================================================

# Function from Adam (but not in AMlib)
fs=sampling
def getOFtemplate(template,noisefft,zeroDC=False,zeroEnds=False):
    '''
    Function to create OF template for Convolution. Includes correct normalization and options to remove DC or zero Ends.
    '''
    tempFFT = np.fft.rfft(template)
    #noisefftN = interNoise(noisefft,template.shape[0])

    noisepsd=noisefft.conjugate()*noisefft
    if zeroDC:
        noisepsd[0]=np.inf
        noisepsd[-1]=np.inf
    tfN = (tempFFT.conjugate() / (noisepsd) ) #/ np.sum(tempFFT.conjugate() * tempFFT / (noisepsd) )
    tfN = np.fft.irfft(tfN)

    #Correct normalization  so divide by max no longer necessary/correct  #tfN = tfN/np.max(tfN)
    zt = tempFFT.conjugate() * tempFFT / noisepsd
    tfN*=tempFFT.shape[-1]/ np.real( np.sum(zt) - zt[0]/2 )  #Subtract zt[0]/2 to correct double counting of f=0 on folded FFT

    if zeroEnds:
        nbins=len(tfN)//20
        tfN+=-np.mean(list(tfN[:nbins])+list(tfN[-nbins:]))
    return tfN

def getBaselineResolution(template,noiseFFT,usePeriodogram=False,fs=fs):
    if usePeriodogram:
        #Skip DC bin [0] which has no bearing on resolution, and is often 0 or inf
        return np.sum(2*scs.periodogram(template,fs=fs)[1][1:] / noiseFFT[1:]**2 )**-0.5  #Periodogram is ^2 , but I always ^(1/2) when handling noiseFFT
    else:
        return np.sum(2*np.abs(np.fft.rfft(template)[1:])**2 / noiseFFT[1:]**2 )**-0.5
    
def psd2fft(psd,fs=fs):
    psdN = psd*1.
    psdN[-1]=psdN[-1]*2**0.5
    return psdN * (fs*(len(psdN)-1))**0.5

#===========================================================

# Define the noise function for shot-noise
def shot_noise(_deltaf,_fb,_pressure,_temperature):

    #bandwidth = np.pi*_fb/2 # Samuli docet
    bandwidth = min(np.pi*_fb/2, lockin_bandwidth) # Samuli docet
    gap = energy_gap_in_low_temp_limit(_pressure)
    mass=mass_effective(_pressure)*atomic_mass # effective mass [kg]
    particle_density=1/(np.pi**2)*np.power(2*mass/Plankbar_const**2,3/2)*np.sqrt(Fermi_momentum(_pressure)**2/(2*mass))*Boltzmann_const*_temperature*np.exp(-gap/(Boltzmann_const*_temperature)) # QP number density, eq.31
    sigma_n=np.sqrt((particle_density*Fermi_velocity(_pressure)*l*diameter/2)/2) # shot-noise, eq.32
    noise = _deltaf/sigma_n * np.sqrt(bandwidth)    # relative error due to the shot-noise in a bandwith???, eq.35

    return noise.item()

#===========================================================


if __name__ == "__main__":

    print('\n==== Hybrid-Toy Simulation ====\n')

    # Parsing arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str, help='Config file', default='config.py')
    parser.add_argument('--noise',type=str, help='Noise model (normal|shot|real)', default='normal')
    args = parser.parse_args()
    exec(open(args.config).read())
    print('Config file:\t',args.config)
    print('Noise model:\t',args.noise)

    # Real settings and noise from Adam's FFT calculated on data for a 'wire'
    if wire==1:
        noiseFFTs = pd.read_csv('data/Bol-1-NoiseFFTs(Volts).csv') # (wrong naming, it contains a sqrt(PSD))
        noiseSum  = pd.read_csv('data/Bol-1-NoiseSummary.csv')
    elif wire==2:
        noiseFFTs = pd.read_csv('data/Bol-1-NoiseFFTs(Volts).csv') # (wrong naming, it contains a sqrt(PSD))
        noiseSum  = pd.read_csv('data/Bol-1-NoiseSummary.csv')

    # Find the row number of the 'time' given as input
    row = noiseSum[noiseSum['Time (s)'] == time].index
    if row.empty:
        print("No row found for the time:", time)
        exit()
    row=row[0] # row number (starting from 0) in the AM/Bol-1-NoiseFFTs(Volts).csv file, corresponding to 2 hours of data
    sumI = noiseSum.iloc[row]
    print('Time:\t',str(int(sumI['Time (s)']))) # time

    # Pressure
    pressure = noiseSum.iloc[row]['Pressure (bar)']
    print('Pressure:\t',pressure,'[bar]')

    # Starting temperature
    print('Starting Temperature:\t',sumI['Temperature (K)']*1e6,'[uK]')
    ttc = sumI['Temperature (K)']/temperature_critical_superfluid(pressure)
    print('Starting T/Tc:\t',ttc)
            
    # Increase of temperature (or T/Tc) vs time
    ttc_rise = ((noiseSum.iloc[row+1]['Temperature (K)']-noiseSum.iloc[row]['Temperature (K)'])/\
            (noiseSum.iloc[row+1]['Time (s)']-noiseSum.iloc[row]['Time (s)']))  # [1/second], increase of T/Tc vs time
    print('Temperature raising:\t',ttc_rise*1e6*3600,'[uK/h]')

    # Real noise
    freqs    = np.array(noiseFFTs['Freqs'])
    noisePSD = np.array(noiseFFTs[str(int(sumI['Time (s)']))]) # sqrt(PSD) [V/sqrt(HZ)] from the Noise.csv file, given a time 
    noisePSD[0]=0
    #noiseFFT=psd2fft(noisePSD) # FFT
    conversion = int(sumI['df/V approx. conversion']) # df/V
    noisePSD = noisePSD*conversion # PSD' [Hz/sqrt[Hz]]
    #noisePSDI = interpolate.interp1d(freqs, noisePSD) # interpolated
    
    # Real parameters for the template
    riseT = sumI['Template Rise Times (1/pi*df) (s)'] # riseT
    fallT = sumI['Template Fall Times (fit) (s)'] + sumI['Template Rise Times (1/pi*df) (s)'] # fallT + riseT
    t_w = riseT
    t_b = fallT + riseT
    
    #if verbose:
    plt.title('Noise PSD from data')
    plt.loglog(freqs,noisePSD)
    plt.xlabel('[Hz]')
    plt.ylabel('Noise PSD $\Delta$f/$\sqrt{Hz}$')
    plt.grid(which='minor',alpha=0.3)
    plt.show()
    
    # Output filename for the true MC
    base, ext = os.path.splitext(filename)
    filename_truth = f"{base}_truth{ext}"
    with open(filename_truth, 'w') as file:
        file.write(f"#time [s] energy [eV]\n")
            
    # full length sample, with correct sampling
    t = np.linspace(0, max_time, max_time*sampling) 
    truth = np.zeros_like(t) # to store truth values

    if args.noise=='real':
        # add noise to the baseline
        total = np.random.normal(0,1,len(t))
        
        original_x = np.linspace(0, 1, len(noisePSD))
        new_x = np.linspace(0, 1, int(len(total)/2+1))
        noisePSD_interpolated = np.interp(new_x, original_x, noisePSD)
        total = np.fft.irfft( np.fft.rfft(total)*psd2fft(noisePSD_interpolated)/len(total)**0.5 )
    else:
        total = np.zeros_like(t)
    
    # random Poisson number of events within the max_time
    N = np.random.poisson(max_time*rate) 
    print('Number of events:\t',N)
    
    # generate N random start times
    starts = np.sort(np.random.randint(0,max_time,N))
        
    # generate a train of N events
    print('\nGenerate events...')
    for start in tqdm(starts):

        # randomised event energy
        energy = np.random.choice(energy_pdf)

        # real increased temperature wrt the starting ttc of the config
        temperature=(ttc+ttc_rise*start)*temperature_critical_superfluid(pressure).item()

        if verbose:
            print('start time [s]',start)
            print('\tenergy',energy)
            print('\ttemperature', temperature)

        # write truth on a file
        with open(filename_truth, 'a') as file:
            file.write(f"{start:.6f} {energy:.9f}\n")
            
        # Base width from the input base temperature
        f_base = Width_from_Temperature(temperature,pressure,diameter)
        
        # Response time
        #t_w = 1/(np.pi*f_base)
        
        delta, _ = DeltaWidth_from_Energy(energy,pressure,temperature,diameter)

        # Winkelmann function: Delta f vs time
        with np.errstate(invalid='ignore'):
            deltaf = np.heaviside(t-start,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-start)/t_b) - np.exp(-(t-start)/t_w)))
        deltaf = np.nan_to_num(deltaf, nan=0.0)  # otherwise there are NaNs before the start of the pulse
        total += deltaf
        truth += deltaf

    #=============================

    # Add non-constant baseline:
    print('\nAdding non-constant baseline...')
    # calculate temperature for each time point
    temperature = (ttc + ttc_rise * t) * temperature_critical_superfluid(pressure)
    # calculate f_base for each temperature
    f_base = np.array([Width_from_Temperature(temp, pressure, diameter) for temp in tqdm(temperature)])
    
    total += f_base  # [Hz]
    truth += f_base  # [Hz]
    
    # Add noise:
    print('\nAdding noise...')
    if args.noise=='normal':
        for i in tqdm(range(len(total))):
            total[i] = np.random.normal(total[i], total[i]/10) # FIX ME
    if args.noise=='shot':
        for i in tqdm(range(len(total))):
            total[i] = np.random.normal(total[i], shot_noise(total[i],f_base[i],pressure,temperature[i]))
        
    # Output:
    print('\nWriting output to',filename,'...')
    with open(filename, 'w') as file:
        file.write(f"#time [s] width [Hz]\n")
        for c1, c2 in zip(t, total):
            file.write(f"{c1:.6f} {c2:.9f}\n")

    # Plotting:
    #if verbose:
    plt.figure(figsize=(15,5))
    plt.plot(t, total, linestyle='',marker='.', color="black", label='Fake data')
    plt.plot(t, truth, linestyle='--', color="green", label='Truth')
    plt.xlabel('time [s]')
    plt.ylabel('$\Delta f$ [Hz]')
    plt.legend()
    plt.show()

    #=======================================

    print('\nAnalysis...')
    
    # Analysis: noise
    
    # template
    tlen = 4096
    taxis = np.arange(tlen)/fs
    template = afp.expon1(taxis,riseT,fallT)
    template = afp.expon1(taxis,riseT,fallT,t_offset=taxis[tlen//2]-taxis[np.argmax(template)])  #Place peak at centre
    template*=1/np.max(template) # normalise
    #template+=-np.mean(template)  # AC Couple

    # Extract noise PSD from the toy trace:
    NTdict = afp.makeNoiseT(total,template,fs=fs,usePeriodogram=True,returnTraces=False,returnPeaks=False, plotHist=False,
                                        #templateThreshold=(10,100), makeTemplateQ=True,tau_rise=riseT,tau_fall=fallT)
                            templateThreshold=(10,100), makeTemplateQ=False,tau_rise=riseT,tau_fall=fallT)
    noisePSD_extracted = NTdict['NoisePSD']
    #noisePSD_extracted = NTdict['NoiseFFT']
    freqs = NTdict['Freqs']

    plt.loglog(freqs,noisePSD, label='Data')
    plt.loglog(freqs,noisePSD_extracted, label='Toy')
    plt.title('Noise PSD comparison')
    plt.xlabel('[Hz]')
    plt.ylabel('Noise PSD $\Delta$f/$\sqrt{Hz}$')
    plt.grid(which='minor',alpha=0.3)
    plt.legend()
    plt.show()

    #=========================================

    # Analysis: peak finder
    
    riseT = sumI['Template Rise Times (1/pi*df) (s)']
    fallT = sumI['Template Fall Times (fit) (s)']
    template = afp.expon1(taxis,riseT,fallT)
    template = afp.expon1(taxis,riseT,fallT,t_offset=taxis[tlen//2]-taxis[np.argmax(template)])  #Place peak at centre
    template*=1/np.max(template)
    if verbose:
        plt.plot(template)
        plt.show()
    
    oftemplate = getOFtemplate(template,psd2fft(noisePSD,fs=fs),zeroDC=True)
    if verbose:
        plt.plot(oftemplate)
        plt.show()
    
    randTraceF = np.convolve(total, oftemplate,mode='valid')
    if verbose:
        plt.plot(randTraceF)
        plt.show()

    # resolution
    calibration = int(sumI['Calibration (Heat+Pulse shape) (keV/Hz)'])  # [MeV/Hz]  # (error in the colum's name in the csv file)
    resB = getBaselineResolution(template,noisePSD_extracted,usePeriodogram=True)  # [Hz]
    print(f'\nBaseline Resolution from toy:  {np.round(resB,3)} [Hz]')
    print(f'Baseline Resolution from toy:  {np.round(resB*calibration*1e6,3)} [eV]')
    res_data=sumI['Measured OF Resolution (eV)']
    print(f'Baseline Resolution from data: {np.round(res_data)} [eV]')

    # peaks
    locs,peak_dict = scs.find_peaks(randTraceF,height = resB*5,prominence=resB*3)
    heights = peak_dict['peak_heights']
    mean,std = np.mean(heights),np.std(heights)
    plt.title(f'Simulated peaks: {N} - Found peaks: {len(locs)}\nExpected resolution: {np.round(resB,2)} [Hz]')
    plt.plot(randTraceF)
    plt.axhline(y=resB*5, color='r', linestyle='--', label='')
    plt.scatter(locs,heights,c='r')
    plt.ylabel('$\Delta f$ [Hz]')
    plt.show()

    plt.hist(heights*calibration, bins=100, log=True)
    plt.xlabel('Energy [MeV]')
    plt.show()

    
