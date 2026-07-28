'''
 - QUEST-DMC WP1 full toy simulation of a train of pulses, for ND3, given
    * a noise FFT    
    * energy PDFs from g4quest and from merged radiogenics
    * energy partition function of measurable events for ER and NR

 - Output:
    Toy simulation single Pickle file:
       * df_total: sample | width with noise [Hz] | true width [Hz] | ID;
       * calibration [eV/Hz];
       * df_truth: ID | start time of the peak [s] | energy [eV] | specie

P. Franchini 7/2026
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

#import matplotlib
#matplotlib.use("Agg")  # use headless backend
#===========================================================

def partition(pdf_file,simulation_events,simulation_rate):
    '''
    Applies the partition of the fraction of measurable energies to a simulation ROOT file
      Arguments: g4quest ROOT file path
      Results: array of energies [keV], merged for ER and NR (WITHOUT the correct total rate normalisation)
    '''

    # Read partition files from E.L. ( Energy [keV] | Measurable fraction | Error )
    partition_ER = pd.read_csv(partition_ER_file)
    partition_ER_energy = partition_ER.iloc[:, 0].to_numpy()
    partition_ER_fraction = partition_ER.iloc[:, 1].to_numpy()
    partition_ER_error = partition_ER.iloc[:, 2].to_numpy()

    partition_NR = pd.read_csv(partition_NR_file)
    partition_NR_energy = partition_NR.iloc[:, 0].to_numpy()
    partition_NR_fraction = partition_NR.iloc[:, 1].to_numpy()
    partition_NR_error = partition_NR.iloc[:, 2].to_numpy()

    # PDG IDs (absolute values)
    pdg_ER = [11, 13, 22]
    pdg_NR = [1000010020, 1000010030, 1000020030, 1000020040]
    pdg_to_be_decided = [211, 2212, 321]  # ???

    pdg_cut_ER = " | ".join(f"(abs(fPDG) == {abs(code)})" for code in pdg_ER)
    pdg_cut_NR = " | ".join(f"(abs(fPDG) == {abs(code)})" for code in pdg_NR)

    # Energy arrays for ER and NR from the ROOT file (geant4 MeV, converted to keV)
    arrays = pdf_file.arrays(["fEdep", "fEvent", "fPDG", "fTrack", "fGlobalTime"], f"(fEdep > 0) & ({pdg_cut_ER})", library="np",)
    pdf_energy_ER = arrays['fEdep']*1e3  # [keV]
    arrays = pdf_file.arrays(["fEdep", "fEvent", "fPDG", "fTrack", "fGlobalTime"], f"(fEdep > 0) & ({pdg_cut_NR})", library="np",)
    pdf_energy_NR = arrays['fEdep']*1e3  # [keV]

    # Interpolate partition fraction vs energy, evaluated on the energy PDFs
    partition_ER_fraction = np.interp(pdf_energy_ER, partition_ER_energy, partition_ER_fraction, left=partition_ER_fraction[0], right=partition_ER_fraction[-1])
    partition_ER_error    = np.interp(pdf_energy_ER, partition_ER_energy, partition_ER_error,    left=partition_ER_error[0], right=partition_ER_error[-1])
    partition_NR_fraction = np.interp(pdf_energy_NR, partition_NR_energy, partition_NR_fraction, left=partition_NR_fraction[0], right=partition_NR_fraction[-1])
    partition_NR_error    = np.interp(pdf_energy_NR, partition_NR_energy, partition_NR_error,    left=partition_NR_error[0], right=partition_NR_error[-1])

    # Random smear generator, keeping the fraction as [0-1]
    rng = np.random.default_rng()
    partition_ER_fraction = np.clip( rng.normal(loc=partition_ER_fraction, scale=partition_ER_error), 0.0, 1.0 )
    partition_NR_fraction = np.clip( rng.normal(loc=partition_NR_fraction, scale=partition_NR_error), 0.0, 1.0 )

    # Partition "convolution" of the energy arrays
    pdf_energy_ER_fraction = pdf_energy_ER * partition_ER_fraction
    pdf_energy_NR_fraction = pdf_energy_NR * partition_NR_fraction

    # Plot comparison before and after the partition is applied in the [0-200keV] range
    if plot:
        # Weights from MC rates
        bins = np.histogram_bin_edges(np.concatenate([pdf_energy_ER]), bins=2000, range=(0,200))
        bin_width = np.diff(bins)[0]  # all bins are uniform, so one is representative
        simulation_weight_per_event = (simulation_rate / simulation_events) * time / bin_width
        simulation_weights = np.full_like(pdf_energy_ER, simulation_weight_per_event)
        if len(pdf_energy_ER>0): plt.hist(pdf_energy_ER, bins=bins, weights=simulation_weights, alpha=0.5, histtype="step",  linewidth=2, label='Simulation', color='orange') # [keV]
        if len(pdf_energy_ER_fraction>0): plt.hist(pdf_energy_ER_fraction, bins=bins, weights=simulation_weights, alpha=0.5, histtype="step",  linestyle='--', linewidth=2, label='Simulation (after partition)') # [keV]
        plt.title('Background simulation - ER')
        plt.xlabel('Deposited energy [keV]')
        plt.yscale('log')
        plt.ylabel('events/day')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        bins = np.histogram_bin_edges(np.concatenate([pdf_energy_NR]), bins=2000, range=(0,200))
        bin_width = np.diff(bins)[0]  # all bins are uniform, so one is representative
        simulation_weight_per_event = (simulation_rate / simulation_events) * time / bin_width
        simulation_weights = np.full_like(pdf_energy_NR, simulation_weight_per_event)
        if len(pdf_energy_NR>0): plt.hist(pdf_energy_NR, bins=bins, weights=simulation_weights, alpha=0.5, histtype="step",  linewidth=2, label='Simulation', color='orange') # [keV]
        if len(pdf_energy_NR_fraction>0): plt.hist(pdf_energy_NR_fraction, bins=bins, weights=simulation_weights, alpha=0.5, histtype="step", linestyle='--', linewidth=2, label='Simulation (after partition)') # [keV]
        plt.title('Background simulation - NR')
        plt.xlabel('Deposited energy [keV]')
        plt.yscale('log')
        plt.ylabel('events/day')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    # Total energies array to be used in the toy
    #print(pdf_energy_ER_fraction)
    #print(pdf_energy_NR_fraction)
    pdf_energy = np.concatenate([pdf_energy_ER_fraction, pdf_energy_NR_fraction])

    return(pdf_energy)


def read_radiogenics(radiogenics): #,simulation_events,simulation_rate):
    '''
    Read a CSV file to generate a PDF for the events injection. Files was produced by the merge code for the radiogenics.
      Arguments: CSV file
      Returns: expected deposited rate [ev/second], (energy_values, energy_probabilities)
    '''
    
    time=86400  # [s] # as in the merged pdf
    
    radiogenics_file = np.loadtxt(radiogenics)
    radiogenics_density = radiogenics_file[:, 0]  # [0-200] keV as in the merged pdf
    radiogenics_energy  = radiogenics_file[:, 1]*1e3  # [ev] as in the merged pdf

    energy_values = radiogenics_energy
    energy_probabilities = radiogenics_density/sum(radiogenics_density) # normalised to 1
    rate = np.sum(radiogenics_density * 0.1)/time

    #HERE add the partition????
    
    return rate, energy_values, energy_probabilities


def read_root(simulation,simulation_events,simulation_rate):
    '''
    Read a g4quest output root file to generate a PDF for the events injection.
    Applies the partition of the fraction of measurable energies.
      Arguments: root file, number of Geant4 (equivalent) primaries, activities [ev/s]
      Returns: expected deposited rate [ev/second], (energy_values [eV], energy_probabilities)
    '''

    # Read Geant4 ROOT file
    pdf_file = uproot.open(f'{simulation}:Scoring')

    # Applies the partition of the fraction of measurable energies
    pdf_energy = partition(pdf_file,simulation_events,simulation_rate)  # array of energies [keV]

    # Weights from MC rates
    bins = np.histogram_bin_edges(np.concatenate([pdf_energy]), bins=10*max_energy, range=(0,max_energy))  # binning of the PDF 0.1 keV/bin
    bin_width = np.diff(bins)[0]  # All bins are uniform, so one is enough
    simulation_weight_per_event = (simulation_rate / simulation_events) * time / bin_width
    simulation_weights = np.full_like(pdf_energy, simulation_weight_per_event)

    hist, bin_edges = np.histogram(pdf_energy, bins=bins, weights=simulation_weights, density=True)
    energy_values = 0.5 * (bin_edges[1:] + bin_edges[:-1]) * 1e3 # [eV]
    energy_probabilities = hist / np.sum(hist)

    rate = len(pdf_energy)*simulation_rate/simulation_events

    return rate, energy_values, energy_probabilities


def inject_events(rate, energy_values, energy_probabilities, truth_ids, description):
    '''
    Inject a train of events in a baseline
    Arguments: rate [ev/second], energy_values, energy_probabilities, truth_ids, description
    Return: truth array
    '''
    # Truth values, full length sample with correct sampling
    t = np.arange(0, max_time, 1/sampling)  # [s]
    truth = np.zeros_like(t) # to store truth values

    # only uses chunks if needed for high max_time (>10 hours)
    if max_time > 36000:
        chunk_size = 36000
    else:
        chunk_size = max_time
    chunks = int(max_time//chunk_size)

    for i in range(chunks):
        chuck_start = i*chunk_size
        chuck_end = chuck_start + chunk_size
        
        # random Poisson number of events within the max_time
        N = np.random.poisson(chunk_size*rate) 
        print('\n'+description)
        print('Number of events to be generated:\t',N)
    
        # generate N random start times, keeping the sampling resolution (1/sampling [s])
        starts = np.sort(np.random.randint(chuck_start*sampling,chuck_end*sampling,N))/100

        # Truth values, full length sample with correct sampling
        t_chunk = np.arange(chuck_start, chuck_end, 1/sampling)  # [s]
        truth_chunk = np.zeros_like(t_chunk) # to store truth values

        # generate a train of N events
        print('Generate events...')
        # Precompute static coefficients
        #exp_factor = np.power(t_b / t_w, t_w / (t_b - t_w))
        #coeff_factor = (t_b / (t_b - t_w))

        for start in tqdm(starts):
         
            # randomised event energy
            energy = np.random.choice(energy_values, p=energy_probabilities)  # generates a random sample from the energy array, based on the 
        
            # real increased temperature wrt the starting ttc of the config
            #temperature=(ttc*start)*temperature_critical_superfluid(pressure).item()

            if verbose:
                print('\ttemperature [K]', temperature)
                print('start time [s]',start)
                print('\tenergy [ev]',energy)
            
            # write truth on a dataframe
            df_truth.loc[len(df_truth)] = [len(df_truth), start, energy, description]  # (ID, start time of the peak [s], energy [eV], species)
        
            # Base width from the input base temperature
            #f_base = Width_from_Temperature(temperature,pressure,diameter)
        
            # Response time
            #t_w = 1/(np.pi*f_base)
            
            #delta, _ = DeltaWidth_from_Energy(energy,pressure,temperature,diameter)
            delta = energy/calibration  # faster for constant temperatures

            if verbose:
                print('\tdelta [Hz]', delta)

            '''
            # Winkelmann function: Delta f vs time (slower method)
            with np.errstate(invalid='ignore'):
            deltaf = np.heaviside(t-start,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-start)/t_b) - np.exp(-(t-start)/t_w)))
            deltaf = np.nan_to_num(deltaf, nan=0.0)  # otherwise there are NaNs before the start of the pulse
            total += deltaf
            truth += deltaf
            '''

            # wire time-constant dependency on the amplitude:
            #t_w = -0.83*delta + 0.18  # from ND3
            
            # Winkelmann function: Delta f vs time (fast method)
            exp_arg1 = -(t_chunk - start) / t_b
            exp_arg2 = -(t_chunk - start) / t_w
            
            # Mask only the valid (t > start) values to avoid NaNs early
            truth_window = 2
            valid_chunk = (t_chunk > start) & (t_chunk <= start + 7*t_b)  # pulse drops < 0.1%
            valid = (t > start) & (t <= start + truth_window*t_b)  # validity for the truth assignment (shorted than the pulse validity to reduce pile up)
            
            deltaf = np.zeros_like(t_chunk)

            # are here not pre-computed
            exp_factor = np.power(t_b / t_w, t_w / (t_b - t_w))
            coeff_factor = (t_b / (t_b - t_w))
            
            coeff = delta * exp_factor * coeff_factor

            deltaf[valid_chunk] = coeff * (np.exp(exp_arg1[valid_chunk]) - np.exp(exp_arg2[valid_chunk]))
            truth_ids[valid] = len(df_truth) - 1  # the index of the truth ID is the 'truth_window', shorter than the pulse

            truth_chunk += deltaf
            
        # Insert chunk results into the main truth array (0, max_time)
        start_idx = int(chuck_start * sampling)
        end_idx = start_idx + len(truth_chunk)
        truth[start_idx:end_idx] = truth_chunk
        
    return truth

def calc_fft_amplitude(t, w):
    '''
    Calculate fft of tracking time series data
    Inputs: time array, width change array
    Outputs: frequency array, fft amplitude array
    E.L.
    '''
    from scipy.fft import fft, fftfreq, rfft, rfftfreq
    import scipy
    
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

def calc_fft_power(t, w):
    '''
    Calculate fft of tracking time series data
    Inputs: time array, width change array
    Outputs: frequency array, fft amplitude array
    E.L.
    '''
    
    t_size = t.size
    s_int = 1/sampling  # sampling interval [s]

    w_noise_fft  = np.fft.rfft(w)
    w_noise_freq = np.fft.rfftfreq(t_size, s_int)
    w_noise_power = (2 / t_size * np.abs(w_noise_fft))**2

    #w_noise_fft = scipy.fftpack.fft(w)
    #w_noise_power = (2 / t_size * np.abs(w_noise_fft))**2
    #w_noise_freq = (scipy.fftpack.fftfreq(t_size, s_int))

    return w_noise_freq, w_noise_power

def shot_noise(_deltaf,_fb,_pressure,_temperature):
    '''
    Define the noise function for shot-noise
    '''
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

    print('\n==== QUEST-DMC Hybrid-Toy Simulation ====\n')

    # Parsing arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str, help='Config file', default='config.py')
    parser.add_argument('--noise',type=str,  help='Noise model (none|normal|shot|real)', default='real')
    parser.add_argument('--output',type=str, help='Output Pickle filename without extension for toy and truth values', default='output')

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
    print('Starting T/Tc:\t\t',ttc)
    f_base = Width_from_Temperature(temperature, pressure, diameter)
    print('Base width [Hz]:\t', f_base)
    
    # Increase of temperature (or T/Tc) vs time
    #ttc_rise = ((noiseSum.iloc[row+1]['Temperature (K)']-noiseSum.iloc[row]['Temperature (K)'])/\
    #        (noiseSum.iloc[row+1]['Time (s)']-noiseSum.iloc[row]['Time (s)']))  # [1/second], increase of T/Tc vs time
    #print('Temperature raising:\t',ttc_rise*1e6*3600,'[uK/h]')

    if args.noise=='real':
        # Load the noise FFT data
        fft_df = pd.read_csv(fft_file)

        #freqs = fft_df['freq [Hz]'].values  # for amplitude
        #amplitudes = fft_df['fft amplitude'].values  # for amplitude
        freqs = fft_df['freq [Hz]'].values[:len(fft_df) // 2]  # for power
        powers = fft_df['power'].values[:len(fft_df) // 2]  # for power
        
        freq_res = freqs[1] - freqs[0]  # frequency resolution from CSV
        t_size_est = int(round(sampling / freq_res))         # estimated t_size
        print('Number of samples of the FFT:\t',t_size_est)

        #fft_rms = np.sqrt(0.5 * np.sum(amplitudes**2))  # for amplitude
        fft_rms = np.sqrt(0.5 * np.sum(np.sqrt(powers)**2))  # for power
        print('Estimated RMS from FFT [mHz]:\t', fft_rms*1e3)

        if plot:
            plt.figure(figsize=(10,4))
            plt.title('Noise FFT from ND3 data - '+fft_file+' - samples: '+str(t_size_est))
            plt.loglog(freqs,powers)
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Power')
            plt.grid(which='minor',alpha=0.3)
            plt.show()
    
    #============================================================

    # Output for the MC truth
    df_truth = pd.DataFrame(columns=['id', 'start', 'energy', 'description'])  # (ID, start time of the peak, energy, species)
    
    t = np.arange(0, max_time, 1/sampling)
    total = np.zeros_like(t)
    truth_ids = np.full_like(t, -1, dtype=int)
    
    calibration = DeltaWidth_from_Energy(1000, pressure, temperature, diameter)[1]  # [eV/Hz]

    print('\nReading ROOT files...')

    cosmics_rate, cosmics_energy_values, cosmics_energy_probabilities = read_root(cosmics,cosmics_events,cosmics_rate)
    source_rate, source_energy_values, source_energy_probabilities = read_root(source,source_events,source_rate)
    gammas_rate, gammas_energy_values, gammas_energy_probabilities = read_root(gammas,gammas_events,gammas_rate)
    neutrons_rate, neutrons_energy_values, neutrons_energy_probabilities = read_root(neutrons,neutrons_events,neutrons_rate)
    if radiogenics_rate>0:
        radiogenics_rate, radiogenics_energy_values, radiogenics_energy_probabilities = read_radiogenics(radiogenics) #,radiogenics_events,radiogenics_rate)

    print('\nInjecting events...')
    if cosmics_rate>0:
        cosmics_truth = inject_events(cosmics_rate, cosmics_energy_values, cosmics_energy_probabilities, truth_ids, 'Cosmics')
        total += cosmics_truth
    if source_rate>0:
        source_truth = inject_events(source_rate, source_energy_values, source_energy_probabilities, truth_ids, 'Source')
        total += source_truth
    if gammas_rate>0:
        gammas_truth = inject_events(gammas_rate, gammas_energy_values, gammas_energy_probabilities, truth_ids, 'Gammas')
        total += gammas_truth
    if neutrons_rate>0:
        neutrons_truth = inject_events(neutrons_rate, neutrons_energy_values, neutrons_energy_probabilities, truth_ids, 'Neutrons')
        total += neutrons_truth
    if radiogenics_rate>0:
        radiogenics_truth = inject_events(radiogenics_rate, radiogenics_energy_values, radiogenics_energy_probabilities, truth_ids, 'Radiogenics')
        total += radiogenics_truth

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    axs[0].hist(
        [df_truth.energy[df_truth.description == 'Cosmics'], df_truth.energy[df_truth.description == 'Source'], df_truth.energy[df_truth.description == 'Gammas'], df_truth.energy[df_truth.description == 'Neutrons'], df_truth.energy[df_truth.description == 'Radiogenics']],
        bins=100,
        label=['Cosmics', 'Source', 'Gammas', 'Neutrons', 'Radiogenics'],
        color=['orange','green','red','blue','darkviolet'],
        histtype='step',
        linewidth=1.5
    )
    axs[0].set_xlim(0, 2e5)
    axs[0].set_yscale('log')
    axs[0].set_title("Simulated energies (full range)")
    axs[0].set_xlabel('Truth Energy [eV]')
    axs[0].legend()
    
    axs[1].hist(
        [df_truth.energy[df_truth.description == 'Cosmics'], df_truth.energy[df_truth.description == 'Source'], df_truth.energy[df_truth.description == 'Gammas'], df_truth.energy[df_truth.description == 'Neutrons'], df_truth.energy[df_truth.description == 'Radiogenics']],
        bins=2000,
        label=['Cosmics', 'Source', 'Gammas', 'Neutrons', 'Radiogenics'],
        color=['orange','green','red','blue','darkviolet'],
        histtype='step',
        linewidth=1.5
    )
    axs[1].set_xlim(0, 1e4)
    axs[1].set_yscale('log')
    axs[1].set_title("Simulated energies (0–10 keV)")
    axs[1].set_xlabel('Truth Energy [eV]')
    
    plt.tight_layout()
    plt.savefig('energy_simulated.png')
    plt.show()
    
    #============================================================

    # Adds a constant baseline:
    total += f_base  # [Hz]
    if cosmics_rate>0: cosmics_truth += f_base  # [Hz]
    if source_rate>0: source_truth += f_base  # [Hz]
    if gammas_rate>0: gammas_truth += f_base  # [Hz]
    if neutrons_rate>0: neutrons_truth += f_base  # [Hz]
    if radiogenics_rate>0: radiogenics_truth += f_base  # [Hz]
    
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
    if args.noise=='none':
        noise = 0
        noisy_trace = total
    if args.noise=='normal':
        noisy_trace = total
        for i in tqdm(range(len(total))):
            #total[i] = np.random.normal(total[i], total[i]/10000)  # FIX ME
            noisy_trace[i] = np.random.normal(total[i], total[i]/10000)  # FIX ME
    if args.noise=='shot':
        for i in tqdm(range(len(total))):
            total[i] = np.random.normal(total[i], shot_noise(total[i],f_base[i],pressure,temperature[i]))
    if args.noise=='real':
        
        '''
        n_samples = int(max_time*sampling)
        
        # 1. compute the frequency bins needed for np.fft.irfft
        #target_fft_freqs = np.fft.rfftfreq(n_samples, d=1/sampling)
        target_fft_freqs = np.fft.rfftfreq(t_size_est, d=1/sampling)

        print("  Original freqs range:   ", freqs[0], freqs[-1])
        print("  Target FFT freqs range: ", target_fft_freqs[0], target_fft_freqs[-1])

        # 2. interpolate amplitude values to match FFT bins
        interp_powers = np.interp(target_fft_freqs, freqs, powers)
        interp_powers = np.sqrt(interp_powers)
        
        # 3. apply random phase to create complex FFT spectrum
        random_phases = np.exp(1j * 2 * np.pi * np.random.rand(len(interp_powers)))
        spectrum = interp_powers * random_phases
        
        # 4. inverse FFT to generate noise in time domain
        noise = np.fft.irfft(spectrum, n=n_samples)
        noise *= n_samples/2
        #noise *= t_size_est/2

        # Scale to the FFT RMS: should be 1, can be removed
        current_rms = np.sqrt(np.mean(noise**2))
        print('  FFT rms scale factor',fft_rms / current_rms)
        noise *= (fft_rms / current_rms)
        '''

        n_samples = int(max_time*sampling)
        #target_fft_freqs = np.fft.rfftfreq(t_size_est, d=1/sampling)
        target_fft_freqs = np.fft.rfftfreq(n_samples, d=1/sampling) # IS THIS A FIX?????
        
        # Original FFT length used to create the CSV
        ##N0 = t.size      # this must be the t_size used in calc_fft_power
        ##fs = sampling

        print("  CSV max freq:       ", freqs[-1])
        print("  Generated max freq: ", target_fft_freqs[-1])
        print("  Lowest CSV freq:", freqs[0])
        print("  Lowest new freq:", target_fft_freqs[1])
        delta_f0 = sampling / t.size

        # ----------------------------------------------------------
        # 1. Convert CSV power to PSD (remove original N dependence)
        psd_values = powers * (t.size / (4 * sampling)) 
        
        # 2. Interpolate PSD to new FFT grid
        interp_psd = np.interp(target_fft_freqs, freqs, psd_values)
        
        # 3. Convert PSD to per-bin power for new n_samples
        delta_f_new = sampling / n_samples
        bin_power = interp_psd * delta_f_new

        # 4. Convert bin power to FFT magnitude
        fft_amplitude = np.sqrt(bin_power) * n_samples
        print("  DC variance contribution:", psd_values[0] * delta_f0)
        print("  Total variance:", np.sum(psd_values * delta_f0))
        
        # Random phase
        random_phases = np.exp(1j * 2 * np.pi * np.random.rand(len(fft_amplitude)))
        spectrum = fft_amplitude * random_phases

        # 5. Inverse FFT
        noise = np.fft.irfft(spectrum, n=n_samples)

        current_rms = np.sqrt(np.mean(noise**2))
        print('  FFT rms scale factor',fft_rms / current_rms)
        noise *= (fft_rms / current_rms)

        '''
        delta_f_new = sampling / n_samples
        variance_generated = np.sum(interp_psd * delta_f_new)
        print("Generated variance:", variance_generated)

        delta_f0 = sampling / t.size
        variance_original = np.sum(psd_values * delta_f0)
        print("Original variance:", variance_original)
        '''
        
        # calculate RMS from the FFT
        power_spectrum = np.abs(spectrum)**2
        power_spectrum[1:-1] *= 2
        noise_fft = np.sqrt(np.sum(power_spectrum) / n_samples**2)
        print('  Estimated RMS from FFT [mHz]:\t', noise_fft*1e3)
        
        # calculate RMS from the generated noise
        noise_rms = np.sqrt(np.mean(noise**2))
        noise_std = np.std(noise)
        print('  Calculated RMS from generated noise [mHz]:\t', noise_rms*1e3)
        print('  Calculated STD from generated noise [mHz]:\t', noise_std*1e3)
        print('  Energy corresponding to the noise RMS [eV]:\t', noise_std*calibration)
        
        plt.figure(figsize=(10, 4))
        plt.title(f'Generated noise from FFT - FFT RMS: {fft_rms:.7f}, Width RMS: {noise_rms:.7f}')
        plt.plot(t, noise)
        plt.xlabel('time [s]')
        plt.ylabel('$\Delta f$ [Hz]')
        plt.show()
        
        # add noise to the signal
        noisy_trace = total + noise
        
        '''
        # Calculate FFT on the noise generated trace
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

        # Compare ND3 noise with the noise generated trace generating PSDs
        # 1. Convert ND3 power spectrum to PSD
        N0 = 2 * len(freqs)          # original FFT length (since freqs = positive half)
        # Convert to PSD [power/Hz]
        psd_csv = powers * (N0 / (4.0 * sampling))

        # 2. Compute PSD of generated noise
        # rFFT of generated noise
        X = np.fft.rfft(noise)
        # True one-sided PSD
        psd_gen = (np.abs(X)**2) / (sampling * noise.size)
        #psd_gen[1:-1] *= 2   # one-sided correction
        freqs_gen = np.fft.rfftfreq(noise.size, d=1/sampling)
        
        plt.figure(figsize=(10,4))
        plt.loglog(freqs, psd_csv, label="PSD from ND3 data")
        plt.loglog(freqs_gen, psd_gen, label="PSD from generated noise")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("PSD [power/Hz]")
        plt.title('Noise FFT comparison - '+fft_file+' - samples: '+str(t_size_est))
        plt.legend()
        plt.grid(True, which="both")
        plt.show()

        '''
        # Calculate FFT on the noise generated trace
        #w_noise_freq, w_noise_amp = calc_fft(t[:t_size_est], noise[:t_size_est]) #t_size_est
        w_noise_freq, w_noise_power = calc_fft_power(t, noise) #t_size_est

        plt.figure(figsize=(10,4))
        plt.title('Noise FFT comparison - '+fft_file+' - samples: '+str(t_size_est))
        plt.loglog(freqs,powers, color='red', alpha=0.7, label='FFT from ND3 data')
        plt.loglog(w_noise_freq,w_noise_power, color='green', alpha=0.7, label='FFT from generated noise')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Power')  # which unit?
        plt.grid(which='minor',alpha=0.3)
        plt.legend()
        plt.show()        
        '''
    
    # Plotting:
    #if verbose:
    plt.figure(figsize=(15,5))
    plt.plot(t, total, linestyle='',marker='.', color='black', label='Fake data')
    plt.plot(t, noisy_trace, linestyle='-', color='red', alpha=0.7, label='Fake data + FFT Noise')
    if cosmics_rate>0:
        plt.plot(t, cosmics_truth, linestyle='--', color='orange', label='Cosmics truth')
    if source_rate>0:
        plt.plot(t, source_truth, linestyle='--', color='green', label='Source truth')
    if gammas_rate>0:
        plt.plot(t, gammas_truth, linestyle='--', color='red', label='Gammas truth')
    if neutrons_rate>0:
        plt.plot(t, neutrons_truth, linestyle='--', color='blue', label='Neutrons truth')
    if radiogenics_rate>0:
        plt.plot(t, radiogenics_truth, linestyle='--', color='blue', label='Radiogenics truth')
    plt.xlabel('time [s]')
    plt.yscale('linear')
    plt.ylabel('$\Delta f$ [Hz]')
    plt.legend(loc='upper right')
    plt.show()

    # Create a DF with `time|width with noise|true width|id`
    df_total = pd.DataFrame({'time': t,'width': noisy_trace,'true_width': total, 'id': truth_ids})  # (t: time, width variation with noise [Hz], true width [Hz], truth_ids: truth ID)
    
    # Output on disk: df_total + calibration + df_truth
    import pickle
    with open(args.output + ".pkl", "wb") as f:
        pickle.dump({"df_total": df_total, "calibration": calibration, "df_truth": df_truth}, f)
    print('\nOutput file:\t',args.output+'.pkl')

    print('\n==== End ====\n')
