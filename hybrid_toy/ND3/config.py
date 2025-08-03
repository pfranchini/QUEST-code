#!/usr/bin/env python

##############################################################
##
##  List of parameters
##
##############################################################

import numpy as np

## Simulation: ###############################################

verbose = False
plot = False

# g4quest output [MeV]
cosmics='/data/questdmc/users/franchinip/QUEST/ND3/cosmics/output-b/cosmics.root'
source='/data/questdmc/users/franchinip/QUEST/ND3/source/output/source.root'

# number of Geant4 (equivalent) primaries
cosmics_events=6.86E+11
source_events=7.58E+11

# rates
cosmics_rate=6e3  # [ev/s], activity*surface of the CRY generator
source_rate=30e3  # [ev/s], 30 kBq Fe55 source

# noise
fft_file='Run23_8mA_01V_noisepwr_quietregions.csv'

max_time = 3600*1   # [second], total lenght of the sample
#max_time = 93040/100   # [second], total lenght of the sample
sampling = 100    # [Hz], sampling (points per second)

filename = "output.txt"  # output's filename

## Cell:  ####################################################

volume = 0.173999e-6  # [m^3] # ND3 bolometer box

## Wire:  ####################################################

diameter = 400e-9  # [m] vibrating wire diameter
l = 2e-3           # [m] vibrating wire length
density = 6.05e3   # [kg/m^3] Niobium-Titanium (NbTi)

t_b = 0.65  # [s] decay constant
t_w = 0.15  # [s] rise constant

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
