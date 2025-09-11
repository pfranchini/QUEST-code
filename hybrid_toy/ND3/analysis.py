'''
 - Simple analysis of the train of pulses:
    * peak finding
    * truth vs reconstruction

Usage example:
    python -i analysis.py --input /data/questdmc/users/franchinip/QUEST/ND3/toy/output

---
P. Franchini 9/2025
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
from scipy.signal import find_peaks
    
# import Tsepelin code
exec(open("../../mod_helium3.py").read())

# import QUEST-DMC code
exec(open("../../mod_quest.py").read())

#===========================================================


if __name__ == "__main__":

    print('\n==== QUEST-DMC Hybrid-Toy Simulation Analysis ====\n')

    # Parsing arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str, help='Config file', default='config.py')
    parser.add_argument('--input',type=str, help='Input Pickle filenames without extension for toy and truth values', default='output')
    
    # Read config file
    args = parser.parse_args()
    exec(open(args.config).read())
    print('Config file:\t',args.config)
    print('Input files:\t',args.input+'*.pkl')

    # TO BE FIXED THE ID OF THE TRUTH THAT IS NOT UNIQUE ACROSS DIFFERENT FILES !!!!
    
    # read Pickle output from toy simulation
    import glob
    files = glob.glob(args.input+'*.pkl')
    print('\nLoading '+str(len(files))+' files...')
    files = [f for f in files if "truth" not in f]
    dfs = [pd.read_pickle(f) for f in files]
    df_total = pd.concat(dfs, ignore_index=True)

    files = glob.glob(args.input+'*_truth.pkl')
    dfs = [pd.read_pickle(f) for f in files]
    df_truth = pd.concat(dfs, ignore_index=True)
    noisy_trace = df_total.width.to_numpy()

    print('\nPeak finding...')
    # compute RMS of noisy trace
    rms_noisy = np.sqrt(np.mean(noisy_trace**2))
    threshold_factor = 1  # define the threshold
    threshold = threshold_factor * rms_noisy  # [Hz]
    
    # find peaks above threshold; peaks are indexes of samples
    peaks, _ = find_peaks(noisy_trace, height=threshold, distance=10*sampling)
    print('Number of peaks: ',len(peaks))

    plt.figure(figsize=(10, 4))
    #plt.plot(total, linestyle='',marker='.', color="black", label='Fake data')
    plt.plot(noisy_trace, label="Fake data + FFT Noise")
    plt.plot(peaks, noisy_trace[peaks], "rx", label=f"Peaks > {threshold_factor}RMS: {len(peaks)}")
    plt.axhline(threshold, color='gray', linestyle='--', label="threshold")
    plt.legend(loc='upper right')
    plt.title("Peak Detection Above RMS Threshold")
    plt.xlabel("Sample")
    plt.ylabel("Width [Hz]")
    plt.savefig('peaks.png')
    plt.show()

    calibration = DeltaWidth_from_Energy(1000, pressure, temperature, diameter)[1]
    energy_threshold = threshold*calibration
    print('Energy threshold [eV]:', energy_threshold)
    print('Min energy detected [eV]:', min(noisy_trace[peaks])*calibration)
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    # Width distribution
    axs[0].hist(noisy_trace[peaks], bins=100)
    axs[0].set_title('Peaks Width Distribution')
    axs[0].set_xlabel('Width [Hz]')
    axs[0].set_yscale('log')
    # Energy distribution
    axs[1].hist(noisy_trace[peaks] * calibration / 1e3, bins=100)
    axs[1].set_title('Energy Distribution')
    axs[1].set_xlabel('Energy [keV]')
    axs[1].set_yscale('log')
    # Energy distribution - zoomed
    axs[2].hist(noisy_trace[peaks] * calibration / 1e3, bins=2000)
    axs[2].set_title('Energy Distribution - zoomed')
    axs[2].set_xlabel('Energy [keV]')
    axs[2].set_xlim(0, 10)
    axs[2].set_yscale('log')
    plt.tight_layout()
    plt.savefig('energy_reconstructed.png')
    plt.show()


    # Reco vs Truth
    # loop through each peak and gather reco and truth energy
    print ('\nReco vs Truth...')
    records = []
    false_positive = 0
    
    for peak in peaks:
        reco_energy = noisy_trace[peak] * calibration
        id_ = df_total.loc[peak, 'id']  # Possible Truth ID on a given peak
        if id_ == -1:
            # it found a fake peak
            false_positive += 1
        truth_row = df_truth[df_truth['id'] == id_]
        if not truth_row.empty:
            truth_energy = truth_row.iloc[0]['energy']
            records.append({'peak': peak, 'id': id_, 'reco_energy': reco_energy, 'truth_energy': truth_energy})

    df_reco_vs_truth = pd.DataFrame(records)
    print(df_reco_vs_truth)
    print('Number of false positives:', false_positive)
    print('Number of false negatives:', len(df_truth) - len(peaks) - false_positive)
    
    plt.figure(figsize=(8,6))
    plt.scatter(df_reco_vs_truth['truth_energy'], df_reco_vs_truth['reco_energy'], alpha=0.7)
    plt.plot([df_reco_vs_truth['truth_energy'].min(), df_reco_vs_truth['truth_energy'].max()],
             [df_reco_vs_truth['truth_energy'].min(), df_reco_vs_truth['truth_energy'].max()], 'r--')
    plt.axvline(x=energy_threshold, color='blue', linestyle='--', linewidth=1.5, label=f'Threshold = {energy_threshold:.1f} eV')
    plt.xlabel("Truth Energy [ev]")
    plt.ylabel("Reconstructed Energy [ev]")
    plt.title("Energy: Reconstructed vs Truth")
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
