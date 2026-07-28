#!/usr/bin/env python

##############################################################
##
##  List of parameters
##
##############################################################

import numpy as np

## Simulation: ###############################################

verbose = False
plot = True

# g4quest output [MeV]
cosmics='/data/questdmc/users/franchinip/QUEST/ND3/cosmics/output-b/cosmics.root'
source='/data/questdmc/users/franchinip/QUEST/ND3/source/output/source.root'
gammas='/home/franchinip/dataQUEST/QUEST/ND3/gammas/output-b/gamma-rethrow.root'  # Fixed Ambient volume
neutrons='/data/questdmc/users/franchinip/QUEST/ND3/neutrons/output-b/neutrons.root'
#neutrons_ambient='/home/franchinip/dataQUEST/QUEST/ND3/neutrons-ambient/output-b/neutrons-ambient.root'
radiogenics='/data/questdmc/users/bloomfield/ND3/QUEST_background_model/output/scaled_summed_0.1kev_rad_pdf.txt'  # Merged; already normalised /cell/day/0.1keV

# number of Geant4 (equivalent) primaries
cosmics_events=6.86E+11
source_events=7.58E+11
gammas_events=1.80E+12  ## old 1.78E+10  # Lizzie
neutrons_events=2.32E+12
neutrons_ambient_events=1.99E+11

# rates
cosmics_rate=6e3  # [ev/s], activity*surface of the CRY generator
source_rate=30e3  # [ev/s], 30 kBq Fe55 source
gammas_rate=2.5*3.66e6 # [ev/s], activity*surface of the ambient generator  ## FIX THIS!!!!
neutrons_rate=6480 # [ev/s], activity*surface of the CRY generator
neutrons_ambient_rate= 0.018*3.66e6 # [ev/s], activity*surface of the ambient generator; activity from Eur. Phys. J. C (2023) 83:94  # FIX THIS
radiogenics_rate=1  # flag to simulate radiogenics (real rate will be read from the file)

# noise
fft_file='Run23_8mA_01V_noisepwr_quietregions.csv'

# energy partition model
partition_NR_file='modelA_NR_response.csv'
partition_ER_file='modelA_NR_response.csv'

time=86400  # [s], time for the binning
max_energy=10000  # [keV], max energy cut

# acquisition
max_time = 3600*10  # [second], total lenght of the sample
#max_time = 93040/100   # [second], total lenght of the sample
sampling = 100  # [Hz], sampling (samples per second)

## Cell:  ####################################################

volume = 0.173999e-6  # [m^3] # ND3 bolometer box

## Wire:  ####################################################

diameter = 400e-9  # [m] vibrating wire diameter
l = 2e-3           # [m] vibrating wire length
density = 6.05e3   # [kg/m^3] Niobium-Titanium (NbTi)

t_b = 0.78  # [s] decay constant from ND3
t_w = 0.15  # [s] rise constant from ND3

pressure = 18.5        # [bar]
temperature = 0.32e-3  # [K]




'''
## Lock-in parameters:  ######################################

amp = 100    # gain of the voltage cold amplifier
v_h = amp*np.pi/2*1e-7  # [V] Base voltage height for a v=1mm/s
lockin_bandwidth = 10   # [Hz] Bandwidth of the lock-in amplifier
v_rms = 7.9*1e-9        # [V] Error on voltage measurement for a lock-in amplifier coming from the datasheet, considering a bandwidth of 10 Hz

## SQUID parameters:  ########################################

B = 0.4e-3  # [T]
R = 1       # [Ohm]
w0 = 5000   # [Hz]
L = 1.5e-6  # [H]
#Z = w0*L  minimum at fields
phi0 = 2.067833848e-15       # [Wb]
S = np.power(0.4e-6*phi0,2)  # [Hz^-1]
M = 10.0e-9  # [H]
v0 = 1e-3    # [m/sec]

##############################################################
'''
