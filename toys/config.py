##############################################################
##
##  List of parameters common to all the toys
##
##############################################################

import numpy as np


## Cell:  ####################################################

volume = 1e-6      # [m^3] Helium-3 cell


## Cryostat:  ################################################

pressure = 0     # [bar] pressure
ttc= 0.1         # T/Tc
energy = 10;     # [eV] deposited energy


## Wire:  ####################################################

diameter = 400e-9;  # [m] vibrating wire diameter
l = 2e-3            # [m] vibrating wire lenght 
density = 6.05e3;   # [kg/m^3] Niobium-Titanium (NbTi)


## Lock-in parameters:  ######################################

t_b = 5.00  # [s] decay constant
amp=100     # gain of the voltage cold amplifier
v_h = amp*np.pi/2*1e-7  # [V] Base voltage height for a v=1mm/s
v_rms = 7.9*1e-9        # [V] Error on voltage measurement for a lock-in amplifier


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
