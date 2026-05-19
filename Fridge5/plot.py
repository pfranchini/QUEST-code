'''
Fridge 5 Bi-207 source simulation
Plot normalised distributions of the deposited energy of the ND3 Geant4 simulations
'''

import numpy as np
import matplotlib.pyplot as plt
import uproot

plt.rcParams.update({'font.size': 14})
#plt.rcParams.update({'lines.linewidth': 3})
plt.rcParams.update({'xtick.direction': 'in'})
plt.rcParams.update({'ytick.direction': 'in'})
plt.rcParams.update({'lines.markersize': 8})

#==========================================================
# g4quest output [MeV]
cosmics='/data/questdmc/users/franchinip/QUEST/ND3/cosmics/output-b/cosmics.root'
neutrons='/data/questdmc/users/franchinip/QUEST/ND3/neutrons/output-b/neutrons.root'
neutrons_ambient='/home/franchinip/dataQUEST/QUEST/ND3/neutrons-ambient/output-b/neutrons-ambient.root'
gammas='/home/franchinip/dataQUEST/QUEST/ND3/gammas/output-b/gamma-rethrow.root'  # Fixed Ambient volume
radiogenics='/data/questdmc/users/bloomfield/ND3/QUEST_background_model/output/scaled_summed_0.1kev_rad_pdf.txt'  # Merged; already normalised /day/0.1keV ?????? - 0-200keV, 2000 bins
source='/data/questdmc/users/franchinip/QUEST/Fridge5/source/output/merge.root'  # Bi-207 simulation for Fridge5

factor=0.230364/0.173999  # rescaling ND3 cell wrt Fridge5 cell volume

# Number of Geant4 (equivalent) primaries
cosmics_events=6.86E+11*factor
neutrons_events=2.32E+12*factor
neutrons_ambient_events=1.99E+11*factor
gammas_events=1.80E+12*factor  ## old 1.78E+10  # Lizzie
source_events=2*100*1E8  # simulation run on Fridge5 (200 files, 1E8 simulated events each)

# Rates
cosmics_rate=6e3 # [ev/s], activity*surface of the CRY generator
neutrons_rate=6480 # [ev/s], activity*surface of the CRY generator
neutrons_ambient_rate= 0.018*3.66e6 # [ev/s], activity*surface of the ambient generator; activity from Eur. Phys. J. C (2023) 83:94  # FIX THIS
gammas_rate=2.5*3.66e6 # [ev/s], activity*surface of the ambient generator  ## FIX THIS!!!!
source_rate=270e3 # [ev/s], 270 kBq Bi207 source in Lancaster

time=86400 # [s]

#==========================================================

cosmics_file = uproot.open(f'{cosmics}:Scoring')
arrays = cosmics_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
cosmics_energy = arrays['fEdep']*1e3  # [keV]

neutrons_file = uproot.open(f'{neutrons}:Scoring')
arrays = neutrons_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
neutrons_energy = arrays['fEdep']*1e3  # [keV]

neutrons_ambient_file = uproot.open(f'{neutrons_ambient}:Scoring')
arrays = neutrons_ambient_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
neutrons_ambient_energy = arrays['fEdep']*1e3  # [keV]

gammas_file = uproot.open(f'{gammas}:Scoring')
arrays = gammas_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
gammas_energy = arrays['fEdep']*1e3  # [keV]

source_file = uproot.open(f'{source}:Scoring')
arrays = source_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
source_energy = arrays['fEdep']*1e3  # [keV]

radiogenics_file = np.loadtxt(radiogenics)

print('Number of deposited events in the ROOT files')
print('  Cosmics         : ', len(cosmics_energy))
print('  Neutrons        : ', len(neutrons_energy))
print('  Neutrons ambient: ', len(neutrons_ambient_energy))
print('  Gammas          : ', len(gammas_energy))
print('  Source          : ', len(source_energy))

print('Expected deposited rate [Hz]')
print('  Cosmics         : ', len(cosmics_energy)*cosmics_rate/cosmics_events)
print('  Neutrons        : ', len(neutrons_energy)*neutrons_rate/neutrons_events)
print('  Neutrons ambient: ', len(neutrons_ambient_energy)*neutrons_ambient_rate/neutrons_ambient_events)
print('  Gammas          : ', len(gammas_energy)*gammas_rate/gammas_events)
print('  Source          : ', len(source_energy)*source_rate/source_events)

print('Expected deposited rate [events/day]')
print('  Cosmics         : ', len(cosmics_energy)*cosmics_rate/cosmics_events*time)
print('  Neutrons        : ', len(neutrons_energy)*neutrons_rate/neutrons_events*time)
print('  Neutrons ambient: ', len(neutrons_ambient_energy)*neutrons_ambient_rate/neutrons_ambient_events*time)
print('  Gammas          : ', len(gammas_energy)*gammas_rate/gammas_events*time)
print('  Source          : ', len(source_energy)*source_rate/source_events*time)

'''
# Plotting ===========

max_energy=200 # [keV]

# Fix the binning for multiple histograms on the same plot
bins = np.histogram_bin_edges(np.concatenate([cosmics_energy, neutrons_energy, source_energy]), bins=100, range=(0,max_energy))
bin_width = np.diff(bins)[0]  # All bins are uniform, so one is enough

# Histogram weights per event
cosmics_weight_per_event = (cosmics_rate / cosmics_events) * time / bin_width
neutrons_weight_per_event = (neutrons_rate / neutrons_events) * time / bin_width
neutrons_ambient_weight_per_event = (neutrons_ambient_rate / neutrons_ambient_events) * time / bin_width
gammas_weight_per_event = (gammas_rate / gammas_events) * time / bin_width
source_weight_per_event   = (source_rate   / source_events)   * time / bin_width

# Weights
cosmics_weights = np.full_like(cosmics_energy, cosmics_weight_per_event)
neutrons_weights = np.full_like(neutrons_energy, neutrons_weight_per_event)
neutrons_ambient_weights = np.full_like(neutrons_ambient_energy, neutrons_ambient_weight_per_event)
gammas_weights = np.full_like(gammas_energy, gammas_weight_per_event)
source_weights = np.full_like(source_energy, source_weight_per_event)

# Plots
plt.hist(cosmics_energy,  bins=bins, weights=cosmics_weights, alpha=0.5, histtype="step",  linewidth=2, label='Cosmics', color='orange')
plt.hist(neutrons_energy, bins=bins, weights=neutrons_weights, alpha=0.5, histtype="step", linewidth=2, label='Neutrons', color='blue')
plt.hist(neutrons_ambient_energy, bins=bins, weights=neutrons_ambient_weights, alpha=0.5, histtype="step", linewidth=2, label='Ambient Neutrons', color='lightblue')
plt.hist(gammas_energy, bins=bins, weights=gammas_weights, alpha=0.5, histtype="step", linewidth=2, label='Gammas', color='red')
plt.hist(source_energy,   bins=bins, weights=source_weights, alpha=0.5, histtype="step", linewidth=2, label='Bi207 Source', color='green')

plt.title('Background simulation')
plt.xlabel('Deposited energy [keV]')
plt.yscale('log')
plt.ylabel('Events/day')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
'''

# Plotting ===========

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

max_energy=1000 # [keV]

bins1 = np.histogram_bin_edges(
    np.concatenate([cosmics_energy, neutrons_energy, neutrons_ambient_energy, gammas_energy, source_energy]), 
    bins=100, 
    range=(0, max_energy)
)
bin_width1 = np.diff(bins1)[0]
print('Bin width',bin_width1)

cosmics_weight_per_event1 = (cosmics_rate / cosmics_events) * time / bin_width1
neutrons_weight_per_event1 = (neutrons_rate / neutrons_events) * time / bin_width1
neutrons_ambient_weight_per_event1 = (neutrons_ambient_rate / neutrons_ambient_events) * time / bin_width1
gammas_weight_per_event1 = (gammas_rate / gammas_events) * time / bin_width1
source_weight_per_event1  = (source_rate / source_events) * time / bin_width1

cosmics_weights1 = np.full_like(cosmics_energy, cosmics_weight_per_event1)
neutrons_weights1 = np.full_like(neutrons_energy, neutrons_weight_per_event1)
neutrons_ambient_weights1 = np.full_like(neutrons_ambient_energy, neutrons_ambient_weight_per_event1)
gammas_weights1 = np.full_like(gammas_energy, gammas_weight_per_event1)
source_weights1 = np.full_like(source_energy, source_weight_per_event1)

radiogenics_density = radiogenics_file[:, 0]  
radiogenics_energy  = radiogenics_file[:, 1]  # [keV]
radiogenics_weights, _ = np.histogram(radiogenics_energy,bins=bins1,weights=radiogenics_density * 0.1)
radiogenics_weights = radiogenics_weights / np.diff(bins1)  # density per unit energy
print("  Radiogenics     : ", np.sum(radiogenics_density * 0.1))

axes[0].hist(cosmics_energy,  bins=bins1, weights=cosmics_weights1, alpha=0.5, histtype="step", linewidth=2, label='Cosmics', color='orange')
axes[0].hist(neutrons_energy, bins=bins1, weights=neutrons_weights1, alpha=0.5, histtype="step", linewidth=2, label='Neutrons', color='blue')
axes[0].hist(neutrons_ambient_energy, bins=bins1, weights=neutrons_ambient_weights1, alpha=0.5, histtype="step", linewidth=2, label='Ambient Neutrons', color='lightblue')
axes[0].hist(gammas_energy,   bins=bins1, weights=gammas_weights1, alpha=0.5, histtype="step", linewidth=2, label='Gammas', color='red')
axes[0].hist(source_energy,   bins=bins1, weights=source_weights1, alpha=0.5, histtype="step", linewidth=2, label='Bi207 Source', color='green')
axes[0].stairs(radiogenics_weights,  bins1, linewidth=2, label='Radiogenics', color='darkviolet')

axes[0].set_title(f'Background Simulation')
axes[0].set_xlabel('Deposited energy [keV]')
axes[0].set_ylabel(f"Events/day/{bin_width1:.0f}keV bin")
axes[0].set_yscale('log')
axes[0].legend()
axes[0].grid(True)

#-----------------------

max_energy = 600  # [keV]
bins2 = np.histogram_bin_edges(
    np.concatenate([cosmics_energy, neutrons_energy, neutrons_ambient_energy, gammas_energy, source_energy]), 
    bins=100, 
    range=(0, max_energy)
)
bin_width2 = np.diff(bins2)[0]
print('Bin width',bin_width2)

cosmics_weight_per_event2 = (cosmics_rate / cosmics_events) * time / bin_width2
neutrons_weight_per_event2 = (neutrons_rate / neutrons_events) * time / bin_width2
neutrons_ambient_weight_per_event2 = (neutrons_ambient_rate / neutrons_ambient_events) * time / bin_width2
gammas_weight_per_event2 = (gammas_rate / gammas_events) * time / bin_width2
source_weight_per_event2  = (source_rate / source_events) * time / bin_width2

cosmics_weights2 = np.full_like(cosmics_energy, cosmics_weight_per_event2)
neutrons_weights2 = np.full_like(neutrons_energy, neutrons_weight_per_event2)
neutrons_ambient_weights2 = np.full_like(neutrons_ambient_energy, neutrons_ambient_weight_per_event2)
gammas_weights2 = np.full_like(gammas_energy, gammas_weight_per_event2)
source_weights2 = np.full_like(source_energy, source_weight_per_event2)

radiogenics_density = radiogenics_file[:100, 0]
radiogenics_energy  = radiogenics_file[:100, 1]  # [keV]
radiogenics_weights, _ = np.histogram(radiogenics_energy,bins=bins2,weights=radiogenics_density * 0.1)
radiogenics_weights = radiogenics_weights / np.diff(bins2)  # density per unit energy

axes[1].hist(cosmics_energy,  bins=bins2, weights=cosmics_weights2, alpha=0.5, histtype="step", linewidth=2, label='Cosmics', color='orange')
axes[1].hist(neutrons_energy, bins=bins2, weights=neutrons_weights2, alpha=0.5, histtype="step", linewidth=2, label='Neutrons', color='blue')
axes[1].hist(neutrons_ambient_energy, bins=bins2, weights=neutrons_ambient_weights2, alpha=0.5, histtype="step", linewidth=2, label='Ambient Neutrons', color='lightblue')
axes[1].hist(gammas_energy,   bins=bins2, weights=gammas_weights2, alpha=0.5, histtype="step", linewidth=2, label='Gammas', color='red')
axes[1].hist(source_energy,   bins=bins2, weights=source_weights2, alpha=0.5, histtype="step", linewidth=2, label='Bi207 Source', color='green')
axes[1].stairs(radiogenics_weights,  bins2,      linewidth=2, label='Radiogenics', color='darkviolet')
axes[0].set_title(f'Background Simulation - Fridge5 (rescaled from ND3) + Bi207')
axes[1].set_xlabel('Deposited energy [keV]')
axes[1].set_ylabel(f"Events/day/{bin_width2:.1f}keV bin")
axes[1].set_yscale('log')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()
