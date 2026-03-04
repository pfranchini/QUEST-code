'''
 - Can read multiple output files of the toy.py code

Usage examples:
    python -i read.py --input /data/questdmc/users/franchinip/QUEST/ND3/toy/output
---
P. Franchini 2/2026
'''

import os
import sys
import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt

#===========================================================


if __name__ == "__main__":

    print('\n==== QUEST-DMC Hybrid-Toy Reader ====\n')

    #===========================================================
    max_time = 3600*100   # [second], total lenght of the sample as in the pkl files
    #===========================================================
    
    # Parsing arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',type=str, help='Input Pickle filenames without extension', default='output')
    
    # Read config file
    args = parser.parse_args()
    #exec(open(args.config).read())
    #print('Config file:\t',args.config)
    print('Input files:\t',args.input+'*.pkl')

    # read Pickle output from toy simulation
    import glob
    files = glob.glob(args.input+'*.pkl')
    print('\nLoading '+str(len(files))+' file(s)...')
    files = [f for f in files if "truth" not in f]
    dfs1 = []
    dfs2 = []
    for file_idx, f in enumerate(tqdm(files), start=1):
        obj = pd.read_pickle(f)
        df1 = obj["df_total"]
        df2 = obj["df_truth"]
        calibration = obj["calibration"]

        df1["id"] = df1["id"].where(df1["id"] == -1, df1["id"] + file_idx * 1_000_000)  # shift IDs by a large offset per file, unless is '-1'
        df2["id"] = df2["id"].where(df2["id"] == -1, df2["id"] + file_idx * 1_000_000)  # shift IDs by a large offset per file, unless is '-1'
        df1["time"] = df1["time"] + max_time*(file_idx-1)  # add to the time shift to order the samples
        df2["start"] = df2["start"] + max_time*(file_idx-1)  # add to the time shift to order the samples; start is not used in the code
        
        dfs1.append(df1)
        dfs2.append(df2)

    print('  Merging DFs...')
    df_total = pd.concat((chunk for chunk in dfs1), ignore_index=True, copy=False)
    df_truth = pd.concat((chunk for chunk in dfs2), ignore_index=True, copy=False)

    noisy_trace = df_total.width.to_numpy()  # width with noise
    print('Total number of injected events:',len(df_truth))

    # add to the df_truth DF the noisy energy
    df_total["energy"] = np.round(df_total["width"]*calibration, 1)
    total = df_total.true_width.to_numpy()  # width without noise

    # plots the simulated energies present in the files
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
    plt.show()

    # plot the width distribution vs time
    sampling=100
    plt.figure(figsize=(16, 4))
    plt.plot(total[:36000*sampling], linestyle='', marker='.', color="black", label='Fake data')
    plt.plot(noisy_trace[:36000*sampling], label="Fake data + FFT Noise")
    plt.legend(loc='upper right')
    plt.title("Width distribution (subset)")
    plt.xlabel("Sample")
    plt.ylabel("Width [Hz]")
    plt.show()
