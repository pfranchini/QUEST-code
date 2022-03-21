'''

Input:
 - Energy

 - Helium-3 Pressure: pressure
 - Wire diameter: d
 - Base temperature: t_base
 - Reponse time; t_w
 - Decay constant: t_b
 - Voltage height: v_h
 - Voltage noise: v_rms

Output:
 - Error on Energy 

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.stats import norm

# import Tsepelin code
execfile("mod_helium3.py")

## Input #####################################################

energy = 1000 # [eV]


## Parameters ################################################

volume = 1e-6      # [m^3] Helium-3 cell
density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)   

#=============================================================

pressure = 10     # [bar] pressure
t_base = 150e-6  # [K] base temperature
d = 200e-9;      # [m] vibrating wire diameter

t_b = 5.0   # [s] decay constant
t_w = 0.77  # [s] response time

v_h=4.6*1e-6   # [V] Base voltage height
v_rms=15*1e-9  # [V] Error on voltage measurement
# v_drive=14.5e-3

N = 1000  # number of toys


## More routines: ###########################################

def Width_from_Temperature(Temperature,PressureBar):
    
    gap = energy_gap_in_low_temp_limit(PressureBar)
    width=np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar)/(2*density*np.pi*d)*np.exp(-gap/(Boltzmann_const*Temperature))
    return width

def Temperature_from_Width(Width,PressureBar):
    
    gap = energy_gap_in_low_temp_limit(PressureBar)
    temperature=-gap/(Boltzmann_const*np.log( Width*2*density*np.pi*d/(np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar)))\
    )
    return temperature

def DeltaWidth_from_Energy(E,PressureBar):
    # Find delta width from the input energy deposition:

    # find fit line for the Width variation vs Deposited energy for the base temperature
    W0=Width_from_Temperature(t_base,PressureBar)
        
    DQ = np.array([])  # delta energy [eV]
    DW = np.array([])  # delta width [Hz]    
    
    for dw in np.arange(0,2.5,0.001):  # Delta(Deltaf)
        T2= Temperature_from_Width(W0+dw,PressureBar)
        T1= Temperature_from_Width(W0,PressureBar)
        DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, PressureBar) - heat_capacity_Cv_B_phase_intergral_from_0(T1, PressureBar)) * volume * 6.242e+18) # [eV]
        DW = np.append(DW,dw)
        
    # Draw the plot 
    plt.plot(DQ/1e3,DW*1e3,label='DQvsDW')
    plt.title('Width variation vs Deposited energy')
    plt.xlim([0, 100])
    plt.xlabel('$\Delta$Q [KeV]')
    plt.ylabel('$\Delta(\Delta f)$ [mHz]')
    plt.show()
    
    # Fit line
    global m
    m, q = np.polyfit(DQ, DW, 1)
    
    # Input delta width from input energy
    deltawidth = E*m
       
    return deltawidth

# Define signal fit function Dw vs time
def df(t,fb,d):
    return fb + np.heaviside(t-5,1)*(d*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-5)/t_b) - np.exp(-(t-5)/t_w)))

##########################################################



if __name__ == "__main__":

    # Input delta(width) from input energy
    delta = DeltaWidth_from_Energy(energy,pressure)
    
    # Base width from the input base temperature
    f_base = Width_from_Temperature(t_base,pressure)

    print("Base width:      ",str(f_base), " Hz")
    print("Width variation: ",delta*1000, " mHz")

    
    t = np.linspace(0, 50, 100) # time

    fit = np.array([]) 

    # Repeat the fit N times
    for i in range(N):
        # Delta f vs time
        deltaf = f_base + np.heaviside(t-5,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-5)/t_b) - np.exp(-(t-5)/t_w)))
        #plt.plot(t, deltaf)
        #plt.show()

        #print("RMS error: ",np.power(deltaf[1],2)/(v_h*f_base)*v_rms)
        # Add noise based on voltage error
        for j in range(len(deltaf)):
            deltaf[j] = np.random.normal(deltaf[j],np.power(deltaf[j],2)/(v_h*f_base)*v_rms, 1)

        # Random noise    
        #noise = np.random.normal(0, 1e-3, len(t))
        #deltaf = deltaf + noise

        #plt.plot(t, deltaf)
        #plt.show()
        
        # Fit the noise'd distribution        
        popt, pcov = curve_fit(df,t,deltaf)

        _, d_fit = popt

        fit = np.append(fit,d_fit/m)
        
    plt.plot(t, deltaf)
    plt.plot(t, df(t,*popt))
    plt.xlabel('time [s]')
    plt.ylabel('$\Delta f$ [Hz]')
    plt.show()

    ## fit
    (mu, sigma) = norm.fit(fit)    
    plt.hist(fit/1000,100)
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Entries')


    
    print(mu,sigma,sigma/mu*100,"%")
    
    
    plt.show()

    
