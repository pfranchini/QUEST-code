#!/usr/bin/env python

##############################################################
##
##  List of parameters
##
##############################################################

import numpy as np

## Simulation: ###############################################

wire = 1 # (1,2)

energy_pdf = np.arange(0, 5.0e6, 1)  # energy spectrum [ev]

rate = 0.01      # [events/second], rate of events
max_time = 36000   # [second], total lenght of the sample
sampling = 100    # [Hz], sampling (points per second)

filename = "output.txt"

verbose = False

## Cell:  ####################################################

volume = 0.315e-6  # [m^3] # Bolometer boxes with 0.5mm hole in 100um PET wall. Estimated time constant 1.0s (https://github.com/slazav/cryo_data/blob/main/f4_data.md)

## Wire:  ####################################################

diameter = 4.5e-6;  # [m] vibrating wire diameter
l = 2e-3            # [m] vibrating wire length
density = 6.05e3;   # [kg/m^3] Niobium-Titanium (NbTi)



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
