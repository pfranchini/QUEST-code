'''
Simulation for the bolometer response with a SQUID readout
considering the resolution as in Lev's note

Input:
 - Energy

 - Helium-3 Pressure: pressure
 - Wire diameter: d
 - Base temperature: t_base
 - central resonating frequency

 - SQUID parameters

Output:
 - Error on Energy 

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.stats import norm

# import Tsepelin code
exec(open("mod_helium3.py").read())


## Parameters ################################################

volume = 1e-6      # [m^3] Helium-3 cell
density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)   

#=============================================================

pressure = 5     # [bar] pressure
t_base = 150e-6  # [K] base temperature
d = 200e-9;      # [m] vibrating wire diameter

# USELESS HERE: t_b = 5.00  # [s] decay constant
#t_w = 0.77  # [s] response time - IS NOT CONSTANT -

# USELESS HERE: v_h = np.pi/2*1e-7  # [V] Base voltage height for a v=1mm/s
# USELESS HERE: v_rms = 3.5*1e-9    # [V] Error on voltage measurement for a lock-in amplifier
# v_drive=14.5e-3

#=============================================================

# USELESS HERE: N = 1000  # number of toys
verbose=False # verbosity for plotting

## SQUID paramenters ========================

B = 0.4e-3 # T
R = 1 # Ohm

w0 = 5000  # Hz
L = 1.5e-6 # H
#Z = w0*L  minimum at fields
phi0 = 2.067833848e-15  # Wb
S = np.power(0.4e-6*phi0,2) # Hz^-1
M = 10.0e-9 # H
l = 1e-3 # m
m = np.power(d/2,2)*np.pi*density # kg/m
v0 = 1e-3 # m/sec

# ===========================================
   

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

def DeltaWidth_from_Energy(E,PressureBar,BaseTemperature):
    # Find delta width from the input energy deposition for a certain base temperature

    # find fit line for the Width variation vs Deposited energy for the base temperature
    W0=Width_from_Temperature(BaseTemperature,PressureBar)
        
    DQ = np.array([])  # delta energy [eV]
    DW = np.array([])  # delta width [Hz]    
    
    #for dw in np.arange(0,2.5,0.001):  # Delta(Deltaf)
    for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)  FASTER
        T2= Temperature_from_Width(W0+dw,PressureBar)
        T1= Temperature_from_Width(W0,PressureBar)
        DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, PressureBar) - heat_capacity_Cv_B_phase_intergral_from_0(T1, PressureBar)) * volume * 6.242e+18) # [eV]
        DW = np.append(DW,dw)
        
    # Draw the plot 
    '''
    if verbose and E==10000: 
        plt.plot(DQ/1e3,DW*1e3,label='DQvsDW')
        plt.title('Width variation vs Deposited energy')
        plt.xlim([0, 100])
        plt.xlabel('$\Delta$Q [KeV]')
        plt.ylabel('$\Delta(\Delta f)$ [mHz]')
        plt.show()
    '''
    
    # Fit line to extract the slope alpha: DQ = alpha * DW
    global alpha
    alpha, _ = np.polyfit(DW, DQ, 1)
    
    # Input delta width from input energy
    deltawidth = E/alpha
       
    return deltawidth, alpha

# SQUID resolution as in Lev's note (dW) [Hz]
def SQUID_Resolution(magnetic_field,resistance): 
    Z=l*np.power(magnetic_field,2)/(2*np.pi*m*f_base)
    V=l*v0*magnetic_field # max voltage across the wire

    return np.sqrt( (np.power(Z+resistance,2) + np.power(w0*L,2))*S*np.pi*f_base/2/np.power(M,2) + 4*Boltzmann_const*t_base*resistance*np.pi*f_base/2 + Boltzmann_const*t_base*l*np.power(magnetic_field,2)/m )/V * f_base
#    return np.sqrt( (np.power(Z+resistance,2) + np.power(w0*L,2))*S*np.pi*f_base/2/np.power(M,2) + 4*Boltzmann_const*t_base*resistance*np.pi*f_base/2 + Boltzmann_const*t_base*l*np.power(magnetic_field,2)/m )/V * f_base 

###########################################################

# error calculation (is not a toy anymore)
def Toy(energy):
    
    print()
    print("Energy:      ",str(energy), " eV")
    # Input delta(width) from input energy
    delta, _ = DeltaWidth_from_Energy(energy,pressure,t_base)
    
    # Base width from the input base temperature
    f_base = Width_from_Temperature(t_base,pressure)

    print("Base width:      ",f_base*1000, " mHz")
    print("Width variation: ",delta*1000,  " mHz")

    # Calculate the resolution from eq.8 in the Lev's note on SQUID
    #for field in np.arange(1e-6,1,1e3): 
    delta_sigma = SQUID_Resolution(B,R)
    
    energy_error=  alpha*delta_sigma


    return energy_error


def Run_Toy(start_energy, end_energy, step):

    global error
    global e

    for energy in np.arange(start_energy, end_energy, step):

        sigma_energy = Toy(energy)
        print("Error:", sigma_energy,"eV", sigma_energy/energy*100,"%")

        error = np.append(error,sigma_energy/energy*100)
        e = np.append(e,energy)

        print(energy,*(sigma_energy/energy*100),file=f)


###########################################################


if __name__ == "__main__":

    # Output file
    f = open("squid_error.txt", "w")
    print("# energy[ev]","error[%]",file=f)

    # Parameters used
    print()
    print("Temperature: ",t_base*1e6, " uk")
    print("Diameter:    ",d*1e9," nm")
    print("Pressure:    ",pressure, "bar")


    # Resolution vs the magnetic field
    energy= 10000
    delta, _ = DeltaWidth_from_Energy(energy,pressure,t_base)
    
    f_base = Width_from_Temperature(t_base,pressure)

    Field = np.array([])  # [T]
    Resolution = np.array([])  # [Hz]
    
    for bfield in np.arange(1e-6, 1, 1e-4):
        Field=np.append(Field,bfield)
        Resolution=np.append(Resolution,SQUID_Resolution(bfield,R))
        
    plt.plot(Field, Resolution)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Magnetic field [T]')
    plt.ylabel('Resolution [Hz]')
    plt.show()

    '''
    # Starts the simulation for a range of energies
   
    global error
    global e
    error = np.array([])
    e = np.array([])

    Run_Toy(100, 900, 50)
    Run_Toy(1000, 9000, 500)
    Run_Toy(10000, 100000, 2000)

    # Plot results
    plt.plot(e/1000,error, linestyle='', marker='o', color="black")
    plt.xscale('log')
    plt.ylim([0, 100])  
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Error [%]')
    plt.savefig('squid_error'+str(d*1e9)+'.pdf')
    plt.savefig('squid_error'+str(d*1e9)+'.png')
    plt.show()


    # Single toy run for a defined energy
    #_ = Toy(10000);

    f.close()

    '''
    
    # Error for differente wire/contact resistance

    for R in (1, 0.1, 0.01, 0.001, 0):

        print("Resistance: ",R)

        error = np.array([])
        e = np.array([])

        Run_Toy(10, 900, 10)
        Run_Toy(1000, 9000, 500)
        Run_Toy(10000, 100000, 2000)

        # Draw a line for each resistance
        plt.plot(e/1000,error,label=str(R)+' $\Omega$')     
    
    # Plot results
    plt.xscale('log')
    plt.ylim([0, 100])  
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Error [%]')
    plt.legend()   
    plt.savefig('squid_error'+str(d*1e9)+'.pdf')
    plt.savefig('squid_error'+str(d*1e9)+'.png')
    plt.show()


    # Single toy run for a defined energy
    #_ = Toy(10000);

    f.close()
        
