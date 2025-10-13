'''
 - Can read multiple output files of the toy code
 - Fit with multiple PDFs both the Truth and the Reco

Usage example:
    python -ifit.py --input /data/questdmc/users/franchinip/QUEST/ND3/toy/output

---
P. Franchini 10/2025
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

from functions import plot_corr_matrix

# import Tsepelin code
exec(open("../../mod_helium3.py").read())

# import QUEST-DMC code
exec(open("../../mod_quest.py").read())

#===========================================================

def fit(energies,pdf1,pdf2):
    
    # --- Define mixture model ---
    def mixture_log_likelihood(params):
        w = params[0]  # mixture weight for PDF1
        if not 0 <= w <= 1:
            return np.inf  # invalid weight
        pdf_vals = w * pdf1(energies) + (1 - w) * pdf2(energies)
        return -np.sum(np.log(pdf_vals + 1e-12))  # negative log-likelihood
    
    # --- Fit the mixture weight ---
    print('  Log likelihood...')
    res = minimize(mixture_log_likelihood, x0=[0.5], bounds=[(0, 1)],method='L-BFGS-B',options={'disp': False})
    w_opt = res.x[0]
    print(f"  PDF1 weight: {w_opt:.3f}")
    print(f"  PDF2 weight: {1 - w_opt:.3f}")
        
    # --- Plots ---
    x = np.linspace(min(energies), max(energies), 500)
    mixture_pdf = w_opt * pdf1(x) + (1 - w_opt) * pdf2(x)
    
    plt.figure(figsize=(8,5))
    plt.hist(energies, bins=100, density=True, alpha=0.4, label='Measured Energy')
    plt.plot(x, w_opt*pdf1(x), 'r--', label='PDF1')
    plt.plot(x, (1-w_opt)*pdf2(x), 'b--', label='PDF2')
    plt.plot(x, mixture_pdf, 'k-', lw=2, label='Mixture fit')
    plt.yscale('log')
    plt.legend()
    plt.xlabel("Energy [eV]")
    plt.ylabel("Probability Density")
    plt.show()

    return(w_opt)

#===========================================================

if __name__ == "__main__":

    print('\n==== QUEST-DMC Hybrid-Toy Simulation Fit ====\n')

    # Parsing arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str, help='Config file', default='config.py')
    parser.add_argument('--input',type=str, help='Input Pickle filenames without extension for toy and truth values', default='output')
    
    # Read config file
    args = parser.parse_args()
    exec(open(args.config).read())
    print('Config file:\t',args.config)
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
    #df_total = pd.concat(dfs1, ignore_index=True, copy=False)
    #df_truth = pd.concat(dfs2, ignore_index=True, copy=False)

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
    print('Total number of true injected events:',len(df_truth))

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


    # Fit TRUTH with PDFs

    print ('\nFitting...')
    from scipy.stats import gaussian_kde
    from scipy.optimize import minimize

    # --- Load energy data ---
    #energies = df_truth["energy"].values
    energies = df_truth[df_truth["energy"]<10000].energy.values

    # --- 1. Load your two PDFs ---
    pdf1_data = pd.read_csv("cosmics_pdf.csv")["Energy [ev]"].values
    pdf2_data = pd.read_csv("source_pdf.csv")["Energy [ev]"].values
    
    ## Estimate continuous PDFs using KDE (Kernel Density Estimation)
    #pdf1 = gaussian_kde(pdf1_data)
    #pdf2 = gaussian_kde(pdf2_data)

    # --- Build normalized histograms for each PDF ---
    bins = np.linspace(min(energies), max(energies), 100)  # same range, fixed bins
    
    hist1, bin_edges = np.histogram(pdf1_data, bins=bins, density=True)
    hist2, _ = np.histogram(pdf2_data, bins=bins, density=True)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    # --- Interpolation functions for lookup ---
    from scipy.interpolate import interp1d
    pdf1 = interp1d(bin_centers, hist1, bounds_error=False, fill_value=0)
    pdf2 = interp1d(bin_centers, hist2, bounds_error=False, fill_value=0)
    
    fit(energies,pdf1,pdf2)
    
    # Peak finding in order to fit reco
    
    print('\nPeak finding...')
    # compute RMS of noisy trace
    rms_noisy = np.sqrt(np.mean(noisy_trace**2))  # [Hz]
    threshold_factor = 1  # define the threshold
    threshold = threshold_factor * rms_noisy  # [Hz]
    threshold = 0.001247408194483884  # [Hz]  TRY THIS from the RMS of the FFT
    
    # find peaks above threshold; peaks are indexes of samples
    #distance=10*t_b*sampling
    from toy import read_root
    #distance=(df_total.iloc[-1].time-df_total.iloc[0].time)/(read_root(source,source_events,source_rate)[0]*max_time*len(files))*sampling  # distance estimate from the rate of injected events given the rate and the real number of simulated events
    #distance=(df_total.iloc[-1].time-df_total.iloc[0].time)/(read_root(gammas,gammas_events,gammas_rate)[0]*max_time*len(files))*sampling  # distance estimate from the rate of injected events given the rate and the real number of simulated events
    #distance=(df_total.iloc[-1].time-df_total.iloc[0].time)/(read_root(cosmics,cosmics_events,cosmics_rate)[0]*max_time*len(files))*sampling  # distance estimate from the rate of injected events given the rate and the real number of simulated events
    distance=(df_total.iloc[-1].time-df_total.iloc[0].time)/(( read_root(cosmics,cosmics_events,cosmics_rate)[0] + read_root(source,source_events,source_rate)[0] )*max_time*len(files))*sampling  # distance estimate from the rate of injected events given the rate and the real number of simulated events
    distance = distance
    print('Distance for the peak finding:', distance,'samples')
    peaks, _ = find_peaks(noisy_trace, height=threshold, distance=distance)
    print('Number of peaks: ',len(peaks))

    if plot:
        print('Plotting...')
        plt.figure(figsize=(16, 4))
        plt.plot(total[:360000*sampling], linestyle='', marker='.', color="black", label='Fake data')
        plt.plot(noisy_trace[:360000*sampling], label="Fake data + FFT Noise")
        plt.plot(peaks[peaks < 360000*sampling], noisy_trace[peaks[peaks < 360000*sampling]], "rx", label=f"Peaks > {threshold_factor}RMS: {len(peaks[peaks < 360000])}")
        #plt.plot(peaks[:360000], noisy_trace[peaks[:360000]], "rx", label=f"Peaks > {threshold_factor}RMS: {len(peaks)}")
        plt.axhline(threshold, color='gray', linestyle='--', label="threshold")
        plt.legend(loc='upper right')
        plt.title("Peak Detection Above RMS Threshold (subset)")
        plt.xlabel("Sample")
        plt.ylabel("Width [Hz]")
        plt.savefig('peaks.png')
        plt.show()
    
    energy_threshold = threshold*calibration
    print('Energy threshold [eV]:', energy_threshold)
    print('Min energy detected [eV]:', min(noisy_trace[peaks])*calibration)

    # Fit peaks one by one with the template
    from scipy.optimize import curve_fit
    
    def template(t, t0, f_base, delta, t_b, t_w):
        H = np.heaviside(t-t0, 1.0)
        factor = delta * np.power(t_b / t_w, t_w / (t_b - t_w)) * (t_b / (t_b - t_w))
        exp_part = np.exp(-(t-t0)/t_b) - np.exp(-(t-t0)/t_w)
        
        return f_base + H * factor * exp_part

    fit_results = []
    time = df_total.time.values
    window= 2  # [s]
    for peak in tqdm(peaks, desc="Processing peaks"):

        i0 = peak - 1*sampling
        i1 = peak + 2*sampling
        xi = time[i0:i1]
        xi_shifted = xi - xi[0]
        yi = noisy_trace[i0:i1]

        # initial guesses
        f00 = time[peak] - (t_b*t_w)/(t_w-t_b)*np.log(t_w/t_b) - xi[0]  # given the peak position is the max
        f_base0 = np.median(yi[:max(1, len(yi)//5)])
        delta0 = np.max(yi) - f_base0
        t_b0, t_w0 = t_b, t_w
        p0 = [f00, f_base0, delta0, t_b0, t_w0]

        try:
            popt, pcov = curve_fit(template, xi_shifted, yi, p0=p0, maxfev=200000)
        except RuntimeError as e:
            print(f"  Curve fitting failed: {e}")
            popt, pcov = None, None
            continue
        except Exception as e:
            print(f"  Unexpected error during curve fitting: {e}")
            popt, pcov = None, None
            continue

        fit_results.append({
            "peak": peak,
            "t0": popt[0],
            "f_base": popt[1],
            "delta": popt[2],
            "t_b": popt[3],
            "t_w": popt[4]
        })
            
        if verbose:
            # plot fitted curve
            yi_fit = template(xi_shifted, *popt)
            yi_test = template(xi_shifted, *p0)
            
            plt.figure(figsize=(6, 3))
            plt.plot(xi, yi, label='Noisy trace', alpha=0.7)
            plt.plot(xi, yi_fit, 'r--', label='Fit')
            plt.xlabel('Time [s]')
            plt.ylabel('Width [Hz]')
            plt.title(f'Peak at index {peak}')
            plt.legend()
            plt.tight_layout()
            plt.show()

    df_peaks = pd.DataFrame(fit_results)

    plot_corr_matrix(pcov)
    
    # Plot t_b vs t_w
    plt.figure()
    plt.scatter(df_peaks["t_b"], df_peaks["t_w"])
    plt.xlabel("t_b")
    plt.ylabel("t_w")
    plt.title("t_b vs t_w")
    plt.show()

    # cuts on the fit results
    peaks = df_peaks[(df_peaks['t_b']<1) & (df_peaks['t_w']<0.20) & (df_peaks['t_b']>0.5) & (df_peaks['t_w']>0.10)].peak

    print('Peaks after fit cuts:', len(peaks))

    df_total = df_total.loc[peaks]
    energies = df_total[df_total["energy"] < 10000].energy.values
    
    fit(energies,pdf1,pdf2)
    
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



    quit()





    
    
    # Plots for True-Positive (TP) and False-Positive (FP) peaks
    df_TP = df_total[df_total["id"] > -1]
    peaks_TP = [p for p in peaks if p in df_TP.index]  # list of TP peaks
    print('Number of TP:', len(peaks_TP))
    
    df_FP = df_total[df_total["id"] == -1]
    peaks_FP = [p for p in peaks if p in df_FP.index]  # list of FP peaks
    print('Number of FP:', len(peaks_FP))

    print('Number of false negatives:', len(df_truth) - len(peaks_TP))
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    # Width distribution
    #axs[0].hist([df_TP.loc[peaks_TP, "width"], df_FP.loc[peaks_FP, "width"], df_truth.energy/calibration],
    axs[0].hist([df_peaks[df_peaks['peak'].isin(peaks_TP)].delta, df_peaks[df_peaks['peak'].isin(peaks_FP)].delta, df_truth.energy/calibration],
                bins=100,
                label=['TP', 'FP', 'truth'],
                color=['green','red','black'],
                histtype='step')
    axs[0].set_title("Peaks Width Distribution")
    axs[0].set_xlabel("Width [Hz]")
    axs[0].set_yscale("log")

    # Energy distribution
    #axs[1].hist([df_TP.loc[peaks_TP, "energy"]/1e3, df_FP.loc[peaks_FP, "energy"]/1e3, df_truth.energy/1e3],
    axs[1].hist([df_peaks[df_peaks['peak'].isin(peaks_TP)].delta*calibration/1e3, df_peaks[df_peaks['peak'].isin(peaks_FP)].delta*calibration/1e3, df_truth.energy/1e3],
                bins=100,
                label=['TP', 'FP', 'truth'],
                color=['green','red','black'],
                histtype='step')
    axs[1].set_title("Energy Distribution")
    axs[1].set_xlabel("Energy [keV]")
    axs[1].set_yscale("log")
    
    # Energy distribution - zoomed
    #axs[2].hist([df_TP.loc[peaks_TP, "energy"]/1e3, df_FP.loc[peaks_FP, "energy"]/1e3, df_truth.energy/1e3],
    axs[2].hist([df_peaks[df_peaks['peak'].isin(peaks_TP)].delta*calibration/1e3, df_peaks[df_peaks['peak'].isin(peaks_FP)].delta*calibration/1e3, df_truth.energy/1e3],
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

    # Save the TP to be used as a PDF
    (df_peaks[df_peaks['peak'].isin(peaks_TP)].delta * calibration).rename('Energy [ev]').to_csv('peaks_calibrated.csv', index=False)
    
    # Reco vs Truth
    # loop through each peak and gather reco and truth energy
    print ('\nReco vs Truth...')
    records = []
    false_positive = 0
    
    for peak in peaks:
        #reco_energy = noisy_trace[peak] * calibration
        reco_energy = df_peaks[df_peaks['peak']==peak].delta.iloc[-1] * calibration  # get the reco from the delta of the template fit
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


    print('Min energy reconstructed [eV]:', np.min(df_reco_vs_truth.reco_energy))
    
    # Energy resolution
    from scipy.stats import norm
    diff = df_reco_vs_truth['reco_energy'] - df_reco_vs_truth['truth_energy']
    mu, std = norm.fit(diff)
    plt.figure(figsize=(7,5))
    n, bins, patches = plt.hist(diff, bins=50, density=True, alpha=0.6, color='skyblue', label='Data')

    # Plot the fitted Gaussian PDF
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 200)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'r-', linewidth=2, label=f'μ={mu:.3f}, σ={std:.3f}')
   
    plt.xlabel("Reco - Truth Energy")
    plt.title("Energy resolution of TP")
    plt.legend()
    plt.grid(True, ls=":")
    plt.show()



    
