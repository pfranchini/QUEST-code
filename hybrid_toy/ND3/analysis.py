'''
 - Can read multiple output files of the toy code
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

    # !!!! FIX THE `time` BECAUSE IS WRONG: NEED TO BE SCALED TO APPEND THE SAMPLES - time is never used in the code, only samples
    # read Pickle output from toy simulation
    import glob
    files = glob.glob(args.input+'*.pkl')
    print('\nLoading '+str(len(files))+' file(s)...')
    files = [f for f in files if "truth" not in f]
    dfs1 = []
    dfs2 = []
    for file_idx, f in enumerate(tqdm(files, desc="  toy files"), start=1):
        obj = pd.read_pickle(f)
        df1 = obj["df_total"]
        df2 = obj["df_truth"]
        calibration = obj["calibration"]
        df1["id"] = df1["id"].where(df1["id"] == -1, df1["id"] + file_idx * 1_000_000)  # shift IDs by a large offset per file, unless is '-1'
        df2["id"] = df2["id"].where(df2["id"] == -1, df2["id"] + file_idx * 1_000_000)  # shift IDs by a large offset per file, unless is '-1'
        dfs1.append(df1)
        dfs2.append(df2)
    df_total = pd.concat(dfs1, ignore_index=True)
    df_truth = pd.concat(dfs2, ignore_index=True)

    '''
    files = glob.glob(args.input+'*_truth.pkl')
    dfs = []
    for file_idx, f in enumerate(tqdm(files, desc="  truth files"), start=1):
        df = pd.read_pickle(f)
        df["id"] = df["id"].where(df["id"] == -1, df["id"] + file_idx * 1_000_000)  # shift IDs by a large offset per file, unless is '-1'
        dfs.append(df)
    df_truth = pd.concat(dfs, ignore_index=True)
    '''
    noisy_trace = df_total.width.to_numpy()  # width with noise
    print('Total number of injected events:',len(df_truth))

    # add to the df_truth DF the noisy energy
    df_total["energy"] = np.round(df_total["width"]*calibration, 1)
    total = df_total.true_width.to_numpy()  # width without noise
    
    # plots the simulated energies present in the files
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    axs[0].hist(
        [df_truth.energy[df_truth.description == 'Cosmics'], df_truth.energy[df_truth.description == 'Source'], df_truth.energy[df_truth.description == 'Gammas']],
        bins=100,
        label=['Cosmics', 'Source', 'Gammas'],
        color=['orange','green','red'],
        histtype='step',
        linewidth=1.5
    )
    axs[0].set_xlim(0, 2e5)
    axs[0].set_yscale('log')
    axs[0].set_title("Simulated energies (full range)")
    axs[0].set_xlabel('Truth Energy [eV]')
    axs[0].legend()

    axs[1].hist(
        [df_truth.energy[df_truth.description == 'Cosmics'], df_truth.energy[df_truth.description == 'Source'], df_truth.energy[df_truth.description == 'Gammas']],
        bins=2000,
        label=['Cosmics', 'Source', 'Gammas'],
        color=['orange','green','red'],
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
    
    print('\nPeak finding...')
    # compute RMS of noisy trace
    rms_noisy = np.sqrt(np.mean(noisy_trace**2))
    threshold_factor = 1  # define the threshold
    threshold = threshold_factor * rms_noisy  # [Hz]
    
    # find peaks above threshold; peaks are indexes of samples
    peaks, _ = find_peaks(noisy_trace, height=threshold, distance=10*sampling)
    print('Number of peaks: ',len(peaks))

    plt.figure(figsize=(10, 4))
    plt.plot(total, linestyle='',marker='.', color="black", label='Fake data')
    plt.plot(noisy_trace, label="Fake data + FFT Noise")
    plt.plot(peaks, noisy_trace[peaks], "rx", label=f"Peaks > {threshold_factor}RMS: {len(peaks)}")
    plt.axhline(threshold, color='gray', linestyle='--', label="threshold")
    plt.legend(loc='upper right')
    plt.title("Peak Detection Above RMS Threshold")
    plt.xlabel("Sample")
    plt.ylabel("Width [Hz]")
    plt.savefig('peaks.png')
    plt.show()

    energy_threshold = threshold*calibration
    print('Energy threshold [eV]:', energy_threshold)
    print('Min energy detected [eV]:', min(noisy_trace[peaks])*calibration)

    '''
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
    '''

    # Plots for True-Positive (TP) and False-Positive (FP) peaks
    df_TP = df_total[df_total["id"] > -1]
    peaks_TP = [p for p in peaks if p in df_TP.index]
    print('Number of TP:', len(peaks_TP))
    
    df_FP = df_total[df_total["id"] == -1]
    peaks_FP = [p for p in peaks if p in df_FP.index]
    print('Number of FP:', len(peaks_FP))

    print('Number of false negatives:', len(df_truth) - len(peaks_TP))
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    # Width distribution
    axs[0].hist([df_TP.loc[peaks_TP, "width"], df_FP.loc[peaks_FP, "width"], df_truth.energy/calibration],
                bins=100,
                label=['TP', 'FP', 'truth'],
                color=['green','red','black'],
                histtype='step')
    axs[0].set_title("Peaks Width Distribution")
    axs[0].set_xlabel("Width [Hz]")
    axs[0].set_yscale("log")

    # Energy distribution
    axs[1].hist([df_TP.loc[peaks_TP, "energy"]/1e3, df_FP.loc[peaks_FP, "energy"]/1e3, df_truth.energy/1e3],
                bins=100,
                label=['TP', 'FP', 'truth'],
                color=['green','red','black'],
                histtype='step')
    axs[1].set_title("Energy Distribution")
    axs[1].set_xlabel("Energy [keV]")
    axs[1].set_yscale("log")
    
    # Energy distribution - zoomed
    axs[2].hist([df_TP.loc[peaks_TP, "energy"]/1e3, df_FP.loc[peaks_FP, "energy"]/1e3, df_truth.energy/1e3],
                bins=2000,
                label=['TP', 'FP', 'truth'],
                color=['green','red','black'],
                histtype='step')
    axs[2].set_title("Energy Distribution - zoomed")
    axs[2].set_xlabel("Energy [keV]")
    axs[2].set_xlim(0, 10)
    axs[2].set_yscale("log")

    plt.legend(loc='upper right')    
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
