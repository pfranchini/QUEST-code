'''
Plot normalised distributions of the deposited energy of the ND3 Geant4 simulations
'''

import numpy as np
import matplotlib.pyplot as plt
import uproot

#==========================================================
# g4quest output [MeV]
cosmics='/data/questdmc/users/franchinip/QUEST/ND3/cosmics/output-b/cosmics.root'
neutrons='/data/questdmc/users/franchinip/QUEST/ND3/neutrons/output-b/neutrons.root'
gammas='/data/questdmc/users/franchinip/QUEST/ND3/gammas/output-b/merge.root'
source='/data/questdmc/users/franchinip/QUEST/ND3/source/output/source.root'

# Number of Geant4 (equivalent) primaries
cosmics_events=6.86E+11
neutrons_events=2.32E+12
gammas_events=1.81E+10 # ????
source_events=7.58E+11

# Rates
cosmics_rate=6e3 # [ev/s], activity*surface of the CRY generator
neutrons_rate=6480 # [ev/s], activity*surface of the CRY generator
gammas_rate=178695 # [ev]s], activity*surface of the ambient generator
source_rate=30e3 # [ev/s], 30 kBq Fe55 source

time=86400 # [s]
max_energy=200 # [keV]

#==========================================================


cosmics_file = uproot.open(f'{cosmics}:Scoring')
arrays = cosmics_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
cosmics_energy = arrays['fEdep']*1e3  # [keV]

neutrons_file = uproot.open(f'{neutrons}:Scoring')
arrays = neutrons_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
neutrons_energy = arrays['fEdep']*1e3  # [keV]

gammas_file = uproot.open(f'{gammas}:Scoring')
arrays = gammas_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
gammas_energy = arrays['fEdep']*1e3  # [keV]

source_file = uproot.open(f'{source}:Scoring')
arrays = source_file.arrays(["fEdep","fEvent", "fPDG","fTrack", "fGlobalTime"], "(fEdep >0)", library = "np")
source_energy = arrays['fEdep']*1e3  # [keV]

print('Number of deposited events')
print('  Cosmics:  ',len(cosmics_energy))
print('  Neutrons: ',len(neutrons_energy))
print('  Gammas:   ',len(gammas_energy))
print('  Source:   ',len(source_energy))

print('Expected deposited rate [Hz]')
print('  Cosmics:  ',len(cosmics_energy)*cosmics_rate/cosmics_events)
print('  Neutrons: ',len(neutrons_energy)*neutrons_rate/neutrons_events)
print('  Gammas:   ',len(gammas_energy)*gammas_rate/gammas_events)
print('  Source:   ',len(source_energy)*source_rate/source_events)

print('Expected deposited rate [events/day]')
print('  Cosmics:  ',len(cosmics_energy)*cosmics_rate/cosmics_events*time)
print('  Neutrons: ',len(neutrons_energy)*neutrons_rate/neutrons_events*time)
print('  Gammas:   ',len(gammas_energy)*gammas_rate/gammas_events*time)
print('  Source:   ',len(source_energy)*source_rate/source_events*time)

# Plotting ===========

# Fix the binning for multiple histograms on the same plot
bins = np.histogram_bin_edges(np.concatenate([cosmics_energy, neutrons_energy, source_energy]), bins=100, range=(0,max_energy))
bin_width = np.diff(bins)[0]  # All bins are uniform, so one is enough

# Histogram weights per event
cosmics_weight_per_event = (cosmics_rate / cosmics_events) * time / bin_width
neutrons_weight_per_event = (neutrons_rate / neutrons_events) * time / bin_width
gammas_weight_per_event = (gammas_rate / gammas_events) * time / bin_width
source_weight_per_event   = (source_rate   / source_events)   * time / bin_width

# Weights
cosmics_weights = np.full_like(cosmics_energy, cosmics_weight_per_event)
neutrons_weights = np.full_like(neutrons_energy, neutrons_weight_per_event)
gammas_weights = np.full_like(gammas_energy, gammas_weight_per_event)
source_weights = np.full_like(source_energy, source_weight_per_event)

# Plots
plt.hist(cosmics_energy,  bins=bins, weights=cosmics_weights, alpha=0.5, histtype="step",  linewidth=2, label='Cosmics', color='orange')
plt.hist(neutrons_energy, bins=bins, weights=neutrons_weights, alpha=0.5, histtype="step", linewidth=2, label='Neutrons', color='blue')
plt.hist(gammas_energy, bins=bins, weights=gammas_weights, alpha=0.5, histtype="step", linewidth=2, label='Gammas', color='red')
plt.hist(source_energy,   bins=bins, weights=source_weights, alpha=0.5, histtype="step", linewidth=2, label='Fe55 Source', color='green')

plt.title('Background simulation')
plt.xlabel('Deposited energy [keV]')
plt.yscale('log')
plt.ylabel('events/day')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# Plotting ===========

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)

max_energy=200 # [keV]

bins1 = np.histogram_bin_edges(
    np.concatenate([cosmics_energy, neutrons_energy, gammas_energy, source_energy]), 
    bins=100, 
    range=(0, max_energy)
)
bin_width1 = np.diff(bins1)[0]
print('Bin width',bin_width1)

cosmics_weight_per_event1 = (cosmics_rate / cosmics_events) * time / bin_width1
neutrons_weight_per_event1 = (neutrons_rate / neutrons_events) * time / bin_width1
gammas_weight_per_event1 = (gammas_rate / gammas_events) * time / bin_width1
source_weight_per_event1   = (source_rate   / source_events)   * time / bin_width1

cosmics_weights1 = np.full_like(cosmics_energy, cosmics_weight_per_event1)
neutrons_weights1 = np.full_like(neutrons_energy, neutrons_weight_per_event1)
gammas_weights1 = np.full_like(gammas_energy, gammas_weight_per_event1)
source_weights1 = np.full_like(source_energy, source_weight_per_event1)

axes[0].hist(cosmics_energy,  bins=bins1, weights=cosmics_weights1, alpha=0.5, histtype="step", linewidth=2, label='Cosmics', color='orange')
axes[0].hist(neutrons_energy, bins=bins1, weights=neutrons_weights1, alpha=0.5, histtype="step", linewidth=2, label='Neutrons', color='blue')
#axes[0].hist(gammas_energy,   bins=bins1, weights=gammas_weights1, alpha=0.5, histtype="step", linewidth=2, label='Gammas', color='red')
axes[0].hist(source_energy,   bins=bins1, weights=source_weights1, alpha=0.5, histtype="step", linewidth=2, label='Fe55 Source', color='green')
axes[0].set_title(f'Background Simulation')
axes[0].set_xlabel('Deposited energy [keV]')
axes[0].set_ylabel('events/day')
axes[0].set_yscale('log')
axes[0].legend()
axes[0].grid(True)

#-----------------------

max_energy = 10  # [keV]
bins2 = np.histogram_bin_edges(
    np.concatenate([cosmics_energy, neutrons_energy, gammas_energy, source_energy]), 
    bins=100, 
    range=(0, max_energy)
)
bin_width2 = np.diff(bins2)[0]
print('Bin width',bin_width2)

cosmics_weight_per_event2 = (cosmics_rate / cosmics_events) * time / bin_width2
neutrons_weight_per_event2 = (neutrons_rate / neutrons_events) * time / bin_width2
gammas_weight_per_event2 = (gammas_rate / gammas_events) * time / bin_width2
source_weight_per_event2   = (source_rate   / source_events)   * time / bin_width2

cosmics_weights2 = np.full_like(cosmics_energy, cosmics_weight_per_event2)
neutrons_weights2 = np.full_like(neutrons_energy, neutrons_weight_per_event2)
gammas_weights2 = np.full_like(gammas_energy, gammas_weight_per_event2)
source_weights2 = np.full_like(source_energy, source_weight_per_event2)

axes[1].hist(cosmics_energy,  bins=bins2, weights=cosmics_weights2, alpha=0.5, histtype="step", linewidth=2, label='Cosmics', color='orange')
axes[1].hist(neutrons_energy, bins=bins2, weights=neutrons_weights2, alpha=0.5, histtype="step", linewidth=2, label='Neutrons', color='blue')
#axes[1].hist(gammas_energy,   bins=bins2, weights=gammas_weights2, alpha=0.5, histtype="step", linewidth=2, label='Gammas', color='red')
axes[1].hist(source_energy,   bins=bins2, weights=source_weights2, alpha=0.5, histtype="step", linewidth=2, label='Fe55 Source', color='green')
axes[0].set_title(f'Background Simulation')
axes[1].set_xlabel('Deposited energy [keV]')
axes[1].set_ylabel('events/day')
axes[1].set_yscale('log')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()
