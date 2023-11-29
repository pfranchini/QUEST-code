#!/usr/bin/env python

##############################################################
##
##  List of parameters common to all the toys
##
##############################################################

import numpy as np


## Cell:  ####################################################

#volume = 1e-6     # [m^3] Helium-3 cell
volume = 0.315e-6  # [m^3] # Bolometer boxes with 0.5mm hole in 100um PET wall. Estimated time constant 1.0s (https://github.com/slazav/cryo_data/blob/main/f4_data.md)


## Cryostat:  ################################################

pressure = 0       # [bar] pressure
ttc = 0.12         # T/Tc
energy = 1000;     # [eV] deposited energy


## Wire:  ####################################################

#diameter = 4.5e-6;  # [m] vibrating wire diameter
diameter = 400e-9;  # [m] vibrating wire diameter
l = 2e-3            # [m] vibrating wire length
density = 6.05e3;   # [kg/m^3] Niobium-Titanium (NbTi)


## Lock-in parameters:  ######################################

t_b = 5.00   # [s] decay constant
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


## Simulation limits:  #######################################

diameter_max = 1000e-9  # [m] vibrating wire max diameter

##############################################################
