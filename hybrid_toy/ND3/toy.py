'''
 - QUEST-DMC WP1 full toy simulation of a train of pulses, for ND3, given
    * a noise FFT    
    * an energy PDF
 - Analysis of the train of pulses:
    * noise
    * peak finding
    * resolution

P. Franchini 7/2025
'''

import os
import sys
import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.signal as scs

# import Tsepelin code
exec(open("../../mod_helium3.py").read())

# import QUEST-DMC code
exec(open("../../mod_quest.py").read())

# import Adam M. library
sys.path.append('../')
import AMlib as afp

#===========================================================

def read_root(simulation,simulation_events,simulation_rate):
    '''
    Read a g4quest output root file to generate a PDF for the events injection.
      Arguments: root file
      Returns: expected deposited rate [ev/second], (energy_values, energy_probabilities)
    '''

    time=86400 # [s]
    max_energy=200 # [keV]
        
    pdf_file = uproot.open(f'{simulation}:Scoring')
    arrays = pdf_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
    pdf_energy = arrays['fEdep']*1e3  # [keV]
    bins = np.histogram_bin_edges(np.concatenate([pdf_energy]), bins=2000, range=(0,max_energy))
    bin_width = np.diff(bins)[0]  # All bins are uniform, so one is enough
    simulation_weight_per_event = (simulation_rate / simulation_events) * time / bin_width
    simulation_weights = np.full_like(pdf_energy, simulation_weight_per_event)    

    if plot:
        plt.hist(pdf_energy,  bins=bins, weights=simulation_weights, alpha=0.5, histtype="step",  linewidth=2, label='Simulation', color='orange') # [keV]
        plt.title('Background simulation')
        plt.xlabel('Deposited energy [keV]')
        plt.yscale('log')
        plt.ylabel('events/day')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    hist, bin_edges = np.histogram(pdf_energy, bins=bins, weights=simulation_weights, density=True)
    energy_values = 0.5 * (bin_edges[1:] + bin_edges[:-1]) * 1e3 # [eV]
    energy_probabilities = hist / np.sum(hist)

    if plot:
        # Plot energy values against weights
        plt.plot(energy_values, energy_probabilities, marker='o', linestyle='-', color='b') # [eV]
        plt.xlabel('Energy [ev]')
        plt.yscale('log')
        plt.ylabel('Entries')
        plt.title('Energy PDF')
        plt.grid(True)
        plt.show()

    rate = len(pdf_energy)*simulation_rate/simulation_events
        
    return rate, energy_values, energy_probabilities


def inject_events(rate, energy_values, energy_probabilities, total, description):
    '''
    Inject a train of events in a baseline
    Arguments: rate, energy_values, energy_probabilities, total
    Return: truth array
    '''
    # random Poisson number of events within the max_time
    N = np.random.poisson(max_time*rate) 
    print('\n'+description)
    print('Number of events to be generated:\t',N)
    
    # generate N random start times
    starts = np.sort(np.random.randint(0,max_time,N))

    # Truth values, full length sample with correct sampling
    t = np.arange(0, max_time, 1/sampling)
    truth = np.zeros_like(t) # to store truth values

    # generate a train of N events
    print('Generate events...')
    for start in tqdm(starts):

        # randomised event energy
        energy = np.random.choice(energy_values, p=energy_probabilities)
    
        # real increased temperature wrt the starting ttc of the config
        #temperature=(ttc*start)*temperature_critical_superfluid(pressure).item()

        if verbose:
            print('start time [s]',start)
            print('\tenergy [ev]',energy)
            print('\ttemperature [K]', temperature)

        # write truth on a file  # REPLACE THIS WITH SOME TABLE TO PLOT AT THE END THE TRUTH 
        with open(filename_truth, 'a') as file:
            file.write(f"{start:.6f} {energy:.9f}\n")
            
        # Base width from the input base temperature
        #f_base = Width_from_Temperature(temperature,pressure,diameter)
        
        # Response time
        #t_w = 1/(np.pi*f_base)
        
        delta, _ = DeltaWidth_from_Energy(energy,pressure,temperature,diameter)

        if verbose:
            print('\tdelta [Hz]', delta)
            
        # Winkelmann function: Delta f vs time
        with np.errstate(invalid='ignore'):
            deltaf = np.heaviside(t-start,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-start)/t_b) - np.exp(-(t-start)/t_w)))
        deltaf = np.nan_to_num(deltaf, nan=0.0)  # otherwise there are NaNs before the start of the pulse
        total += deltaf
        truth += deltaf

    return truth

def calc_fft(t, w):
    '''
    Calculate fft of tracking time series data
    Inputs: time array, width change array
    Outputs: frequency array, fft amplitude array
    E.L.
    '''
    from scipy.fft import fft, fftfreq, rfft, rfftfreq
    
    t_size = t.size
    s_int = 1/sampling  # sampling interval, s

    #w_noise_fft = scipy.fftpack.fft(w)
    w_noise_fft = fft(w)
    w_noise_fft = rfft(w) # only real, avoid negative frequencies
    w_noise_amp = 2 / t_size * np.abs(w_noise_fft)
    #w_noise_freq = np.abs(scipy.fftpack.fftfreq(t_size, s_int))
    w_noise_freq = fftfreq(t_size, d=s_int)
    w_noise_freq = rfftfreq(t_size, d=s_int) # only real
   
    return w_noise_freq, w_noise_amp

#----

def getBaselineResolution(template,noiseFFT,usePeriodogram,fs):
    if usePeriodogram:
        #Skip DC bin [0] which has no bearing on resolution, and is often 0 or inf
        return np.sum(2*scs.periodogram(template,fs=fs)[1][1:] / noiseFFT[1:]**2 )**-0.5  #Periodogram is ^2 , but I always ^(1/2) when handling noiseFFT
    else:
        return np.sum(2*np.abs(np.fft.rfft(template)[1:])**2 / noiseFFT[1:]**2 )**-0.5
    
def psd2fft(psd,fs):
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
    parser.add_argument('--noise',type=str, help='Noise model (none|normal|shot|real)', default='normal')

    # Read config file
    args = parser.parse_args()
    exec(open(args.config).read())
    print('Config file:\t',args.config)
    print('Noise model:\t',args.noise)
    
    exec(open(args.config).read())
    print('Diameter:\t',diameter*1e9,'[nm]')
    print('Pressure:\t',pressure,'[bar]')

    # Starting temperature
    print('Starting Temperature:\t',temperature*1e6,'[uK]')
    ttc = temperature/temperature_critical_superfluid(pressure)
    print('Starting T/Tc:\t',ttc)
    f_base = Width_from_Temperature(temperature, pressure, diameter)
    print('Base width [Hz]:\t', f_base)
    
    # Increase of temperature (or T/Tc) vs time
    #ttc_rise = ((noiseSum.iloc[row+1]['Temperature (K)']-noiseSum.iloc[row]['Temperature (K)'])/\
    #        (noiseSum.iloc[row+1]['Time (s)']-noiseSum.iloc[row]['Time (s)']))  # [1/second], increase of T/Tc vs time
    #print('Temperature raising:\t',ttc_rise*1e6*3600,'[uK/h]')

    if args.noise=='real':
        # Load the noise FFT data
        fft_df = pd.read_csv(fft_file)
        freqs = fft_df['freq [Hz]'].values
        amplitudes = fft_df['fft amplitude'].values

        freq_res = freqs[1] - freqs[0]  # frequency resolution from CSV
        # estimated t_size
        t_size_est = int(round(sampling / freq_res))
        print('Number of samples of the FFT:\t',t_size_est)

        fft_rms = np.sqrt(0.5 * np.sum(amplitudes**2))
        print('Estimated RMS from FFT:\t', fft_rms)

        if plot:
            plt.figure(figsize=(10,4))
            plt.title('Noise FFT from ND3 data - '+fft_file+' - samples: '+str(t_size_est))
            plt.loglog(freqs,amplitudes)
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Amplitude')
            plt.grid(which='minor',alpha=0.3)
            plt.show()
    
    # Output filename for the true MC
    base, ext = os.path.splitext(filename)
    filename_truth = f"{base}_truth{ext}"
    with open(filename_truth, 'w') as file:
        file.write(f"#time [s] energy [eV]\n")

    #t = np.linspace(0, max_time, max_time*sampling)
    t = np.arange(0, max_time, 1/sampling)
    total = np.zeros_like(t)

    #============================================================
        
    cosmics_rate, cosmics_energy_values, cosmics_energy_probabilities = read_root(cosmics,cosmics_events,cosmics_rate)
    source_rate, source_energy_values, source_energy_probabilities = read_root(source,source_events,source_rate)

    energy_values=cosmics_energy_values
    energy_probabilities=cosmics_energy_probabilities

    cosmics_truth = inject_events(cosmics_rate, cosmics_energy_values, cosmics_energy_probabilities, total, 'Cosmics')
    source_truth = inject_events(source_rate, source_energy_values, source_energy_probabilities, total, 'Ion-55 source')
    
    #=============================

    # Adds a constant baseline:
    #total += f_base  # [Hz]
    #cosmics_truth += f_base  # [Hz]
    #source_truth += f_base  # [Hz]

    '''
    # Add non-constant baseline:    # Skipping for now since not much effect in 2 hour traces
    print('\nAdding non-constant baseline...')
    # calculate temperature for each time point
    temperature = (ttc + ttc_rise * t) * temperature_critical_superfluid(pressure)
    # calculate f_base for each temperature
    f_base = np.array([Width_from_Temperature(temp, pressure, diameter) for temp in tqdm(temperature)])
    total += f_base  # [Hz]
    truth += f_base  # [Hz]
    '''
    # Add noise:
    print('\nAdding noise...')
    if args.noise=='normal':
        for i in tqdm(range(len(total))):
            total[i] = np.random.normal(total[i], total[i]/10000)  # FIX ME
    if args.noise=='shot':
        for i in tqdm(range(len(total))):
            total[i] = np.random.normal(total[i], shot_noise(total[i],f_base[i],pressure,temperature[i]))
    if args.noise=='real':

        n_samples = int(max_time*sampling)
        
        # 1. compute the frequency bins needed for np.fft.irfft
        #target_fft_freqs = np.fft.rfftfreq(n_samples, d=1/sampling)
        target_fft_freqs = np.fft.rfftfreq(t_size_est, d=1/sampling)

        #print("Original freqs range: ", freqs[0], freqs[-1])
        #print("Target FFT freqs range: ", target_fft_freqs[0], target_fft_freqs[-1])

        # 2. interpolate amplitude values to match FFT bins
        interp_amplitudes = np.interp(target_fft_freqs, freqs, amplitudes)
        #interp_amplitudes = amplitudes
        
        # 3. apply random phase to create complex FFT spectrum
        random_phases = np.exp(1j * 2 * np.pi * np.random.rand(len(interp_amplitudes)))
        spectrum = interp_amplitudes * random_phases
        
        # 4. inverse FFT to generate noise in time domain
        noise = np.fft.irfft(spectrum, n=n_samples)
        noise *= n_samples/2
        #noise *= t_size_est/2
        
        #current_rms = np.sqrt(np.mean(noise**2))
        #fft_rms = 0.0009700914070394942 # quite RMS from EL
        #noise *= (fft_rms / current_rms)  # Scale to the FFT RMS

        # calculate RMS from the FFT
        power_spectrum = np.abs(spectrum)**2
        power_spectrum[1:-1] *= 2
        noise_fft = np.sqrt(np.sum(power_spectrum) / n_samples**2)

        # calculate RMS from the generated noise
        noise_rms = np.sqrt(np.mean(noise**2))

        plt.figure(figsize=(10, 4))
        plt.title(f'Generated noise from FFT - FFT RMS: {fft_rms:.5f}, Width RMS: {noise_rms:.5f}')
        plt.plot(t, noise)
        plt.xlabel('time [s]')
        plt.ylabel('$\Delta f$ [Hz]')
        plt.show()
        
        # add noise to the signal
        noisy_trace = total + noise

        '''
        # Compute FFT of generated noise
        noise_fft = np.fft.rfft(noise)
        noise_amplitude_spectrum = np.abs(noise_fft)
        noise_freqs = np.fft.rfftfreq(len(noise), d=1/sampling)

        plt.figure(figsize=(12, 5))
        plt.loglog(freqs, amplitudes, label='Original Spectrum', linewidth=2)
        plt.loglog(noise_freqs, noise_amplitude_spectrum, label='FFT of Generated Noise', alpha=0.7)
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Amplitude")
        plt.title("Comparison: Original Spectrum vs. FFT of Generated Noise")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        '''

        # Calculate FFT on the noise generated trace
        w_noise_freq, w_noise_amp = calc_fft(t[:t_size_est], noise[:t_size_est]) #t_size_est
        w_noise_freq, w_noise_amp = calc_fft(t, noise) #t_size_est
        
        plt.figure(figsize=(10,4))
        plt.title('Noise FFT comparison - '+fft_file)
        plt.loglog(freqs,amplitudes, alpha=0.7, label='FFT from ND3 data')
        plt.loglog(w_noise_freq, w_noise_amp, alpha=0.7, label='FFT from generated noise')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude')
        plt.grid(which='minor',alpha=0.3)
        plt.legend()
        plt.show()        

    
    # Plotting:
    #if verbose:
    plt.figure(figsize=(15,5))
    plt.plot(t, total, linestyle='',marker='.', color="black", label='Fake data')
    plt.plot(t, noisy_trace, linestyle='-', color="red", alpha=0.7, label='Fake data + FFT Noise')
    plt.plot(t, cosmics_truth, linestyle='--', color="orange", label='Cosmics truth')
    plt.plot(t, source_truth, linestyle='--', color="green", label='Source truth')
    plt.xlabel('time [s]')
    plt.yscale('linear')
    plt.ylabel('$\Delta f$ [Hz]')
    plt.legend()
    plt.show()

    # Output:
    print('\nWriting output to',filename,'...')
    with open(filename, 'w') as file:
        file.write(f"#time [s] width [Hz]\n")
        for c1, c2 in zip(t, total):
            file.write(f"{c1:.6f} {c2:.9f}\n")

    #================================================

    print('\nAnalysis...')

    import numpy as np
    from scipy.signal import find_peaks
    import matplotlib.pyplot as plt

    # 1. Compute RMS of noisy trace
    rms_noisy = np.sqrt(np.mean(noisy_trace**2))
    
    # 2. Define threshold factor (e.g., 3×RMS)
    threshold_factor = 3
    threshold = threshold_factor * rms_noisy
    
    # 3. Find peaks above threshold
    peaks, _ = find_peaks(noisy_trace, height=threshold, distance=10*sampling)

    print('Number of peaks: ',len(peaks))
        
    # 4. Plot
    plt.figure(figsize=(10, 4))
    plt.plot(noisy_trace, label="Noisy Trace")
    plt.plot(peaks, noisy_trace[peaks], "rx", label=f"Peaks > {threshold_factor}RMS: {len(peaks)}")
    plt.axhline(threshold, color='gray', linestyle='--', label="threshold")
    plt.legend()
    plt.title("Peak Detection Above RMS Threshold")
    plt.xlabel("Sample")
    plt.ylabel("Width [Hz]")
    plt.show()

    plt.hist(noisy_trace[peaks], bins=100, log=False)
    plt.title('Peaks width distribution')
    plt.xlabel('Width [Hz]')
    plt.show()
    
    
    quit()



    



    
    # Old analysys of F#4 to be removed...
    
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
    resB = getBaselineResolution(template,noisePSD_extracted,usePeriodogram=True,fs=fs)  # [Hz]
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

    
