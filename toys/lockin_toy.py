'''
Toy MC simulation for the bolometer response
using the error from the LOCKIN circuit

Input:

 - Energy
 - Helium-3 Pressure: pressure
 - Wire diameter: diameter
 - Base temperature: temperature

 - Reponse time; t_w
 - Decay constant: t_b
 - Voltage height: v_h
 - Voltage noise: v_rms

Output:
 - Preliminary plots: voltage error
 - Error vs Pressure, Temperature, Diameter

P. Franchini, 10/2022

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.stats import norm

# import Tsepelin code
exec(open("../mod_helium3.py").read())

## Parameters ################################################

volume = 1e-6      # [m^3] Helium-3 cell
density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)   

#=============================================================

pressure = 5     # [bar] pressure
#temperature = 150e-6    # [K] base temperature - useless, we need just T/Tc
ttc=0.1  # T/Tc
diameter = 200e-9;      # [m] vibrating wire diameter
energy = 10;   # [eV] deposited energy

#=============================================================

t_b = 5.00  # [s] decay constant
amp=100 # gain of the voltage cold amplifier
v_h = amp*np.pi/2*1e-7  # [V] Base voltage height for a v=1mm/s
v_rms = 7.9*1e-9    # [V] Error on voltage measurement for a lock-in amplifier

#=============================================================

N = 1000  # number of toys
verbose=False  # verbosity for plotting

unused = 0.0
## LOCKIN parameters ==========================================                                                                                                            

# ===========================================================

## More routines: ###########################################

def Width_from_Temperature(Temperature,PressureBar,Diameter):
    
    gap = energy_gap_in_low_temp_limit(PressureBar)
    width=np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar)/(2*density*np.pi*Diameter)*np.exp(-gap/(Boltzmann_const*Temperature))
    return width

def Temperature_from_Width(Width,PressureBar,Diameter):
    
    gap = energy_gap_in_low_temp_limit(PressureBar)
    temperature=-gap/(Boltzmann_const*np.log( Width*2*density*np.pi*Diameter/(np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar))))
    return temperature

def DeltaWidth_from_Energy(E,PressureBar,BaseTemperature,Diameter):
    # Find delta width from the input energy deposition for a certain base temperature

    # find fit line for the Width variation vs Deposited energy for the base temperature
    W0=Width_from_Temperature(BaseTemperature,PressureBar,Diameter)

    DQ = np.array([])  # delta energy [eV]
    DW = np.array([])  # delta width [Hz]    
    
    #for dw in np.arange(0,2.5,0.001):  # Delta(Deltaf)
    for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)  FASTER
        T2= Temperature_from_Width(W0+dw,PressureBar,Diameter)
        T1= Temperature_from_Width(W0,PressureBar,Diameter)
        DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, PressureBar) - heat_capacity_Cv_B_phase_intergral_from_0(T1, PressureBar)) * volume * 6.242e+18) # [eV]
        DW = np.append(DW,dw)
        
    # Fit line to extract the slope alpha: DQ = alpha * DW
    global alpha
    alpha, _ = np.polyfit(DW, DQ, 1)
    
    # Input delta width from input energy
    deltawidth = E/alpha
       
    return deltawidth, alpha


# Define signal fit function Dw vs time
def df(_t,_fb,_d): # time, base width, delta (delta width)
    _t_w = 1/(np.pi*_fb)
    _t1=5
    return _fb + np.heaviside(_t-_t1,1)*(_d*np.power(t_b/_t_w,_t_w/(t_b-_t_w))*(t_b/(t_b-_t_w))*(np.exp(-(_t-_t1)/t_b) - np.exp(-(_t-_t1)/_t_w)))

###########################################################

def Toy(_energy,_pressure,_temperature,_diameter):

    print()
    # Input delta(width) from input energy
    delta, _ = DeltaWidth_from_Energy(_energy,_pressure,_temperature,_diameter)
    
    # Base width from the input base temperature
    f_base = Width_from_Temperature(_temperature,_pressure,_diameter)

    # Response time
    t_w = 1/(np.pi*f_base)

    print("Pressure:    ",_pressure, " bar")
    print("Temperature: ",_temperature*1e6, " uk")
    print("Diameter:    ",_diameter*1e9, " nm")
    print("T/Tc:        ",_temperature/temperature_critical_superfluid(_pressure))

    print("Base width:      ",f_base*1000, " mHz")
    print("Width variation: ",delta*1000,  " mHz")
    print("t_w: ",t_w, "s")
    print("Bandwidth: ",np.pi*f_base/2, " Hz")
    
    t = np.linspace(4.5, 50, 1000) # time

    base_toy = np.array([])  # base width distribution
    delta_toy = np.array([]) # delta width distribution

    # Repeat the fit N times
    for i in range(N):

        # Delta f vs time
        #print(np.exp(-(t-5)/t_b))
        #print(np.exp(-(t-5)/t_w))
        deltaf = f_base + np.heaviside(t-5,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-5)/t_b) - np.exp(-(t-5)/t_w)))
        
        # Add noise based on voltage error. The voltage error comes from the LOCKIN circuit.
        for j in range(len(deltaf)):
            deltaf[j] = np.random.normal(deltaf[j],np.power( deltaf[j],2 ) / (v_h*f_base)*v_rms, 1)

        # Fit the noise'd distribution        
        popt, pcov = curve_fit(df,t,deltaf)
        # errors on base width and width increase
        base_fit, delta_fit = popt

        delta_toy = np.append(delta_toy,delta_fit)
        base_toy = np.append(base_toy,base_fit)

        
    if (verbose):
        # Plot deltaf(t)                                                                                                                                                   
        plt.title(str(_pressure)+' bar - '+str(diameter*1e9)+' nm - '+str(_temperature*1e6)+" $\mu$K")
        plt.plot(t, deltaf,linestyle='',marker='.', color="black")
        plt.plot(t, df(t,*popt))
#        plt.ylim([0,delta*1.1])
        plt.xlabel('time [s]')
        plt.ylabel('$\Delta f$ [Hz]')
        plt.show()
    
    # Plot toy energy distribution
    if verbose and _energy==10000:
        plt.hist(delta_toy*alpha/1000,100)
        plt.xlabel('Energy [KeV]')
        plt.ylabel('Entries')
        plt.savefig('energy_distribution'+str(d*1e9)+'.pdf')
        plt.show()

    # Plot toy base distribution
    if verbose and _energy==10000:
        plt.hist(base_toy,100)
        plt.xlabel('Base width [Hz]')
        plt.ylabel('Entries')
        plt.savefig('base_width_distribution'+str(d*1e9)+'.pdf')
        plt.show()


    ## Gaussian fit for base and delta toy distributions
    (delta_mu, delta_sigma) = norm.fit(delta_toy)
    (base_mu, base_sigma) = norm.fit(base_toy)

    print("Fitted base: ",base_mu," ",base_sigma)
    print("Fitted delta: ",delta_mu," ",delta_sigma)

    # Calculates alpha_prime for the error propagation
    epsilon=1e-9
    _, alpha1 = DeltaWidth_from_Energy(_energy, _pressure, _temperature - epsilon/2, _diameter)
    _, alpha2 = DeltaWidth_from_Energy(_energy, _pressure, _temperature + epsilon/2, _diameter)
    gap = energy_gap_in_low_temp_limit(_pressure)
    alpha_prime = (alpha2-alpha1)/epsilon * Boltzmann_const/gap * np.power(_temperature,2) * 1/Width_from_Temperature(_temperature,_pressure,_diameter)
    print("alpha_prime",alpha_prime)
  
    print(alpha_prime*delta_mu*base_sigma)
    print(alpha*delta_sigma)
    
    _error=  np.sqrt( np.power(alpha_prime*delta_mu*base_sigma,2) + np.power(alpha*delta_sigma,2) )
    #energy_error=  np.sqrt( np.power(alpha_prime*delta_mu*base_sigma,2) )

    return _error


def Run_Toy_Pressure(start, end, step, _pressure,_ttc,_diameter, _f):

    global error
    global value

    for v in np.arange(start, end, step):
        _temperature=_ttc*temperature_critical_superfluid(v)
        sigma = Toy(energy,v,_temperature,_diameter)
        print("Error:", sigma,"eV, ", sigma/energy*100,"%")

        error = np.append(error,sigma/energy*100)
        value = np.append(value,v)

        print(v,*(sigma/energy*100),file=_f)

        
def Run_Toy_Diameter(start, end, step, _pressure,_ttc,_diameter, _f):

    global error
    global value

    _temperature=_ttc*temperature_critical_superfluid(_pressure)
    for v in np.arange(start, end, step):
        sigma = Toy(energy,_pressure,_temperature,v)
        print("Error:", sigma,"eV, ", sigma/energy*100,"%")

        error = np.append(error,sigma/energy*100)
        value = np.append(value,v)

        print(v,*(sigma/energy*100),file=_f)


def Run_Toy_Temperature(start, end, step, _pressure,_temperature,_diameter, _f):

    global error
    global value

    for v in np.arange(start, end, step):
        sigma = Toy(energy,_pressure,v*temperature_critical_superfluid(_pressure),_diameter)
        
        print("Error:", sigma,"eV, ", sigma/energy*100,"%")

        error = np.append(error,sigma/energy*100)
        value = np.append(value,v)

        print(v,*(sigma/energy*100),file=_f)

                
###########################################################


if __name__ == "__main__":

    # Output files
    f1 = open("output/lockin_toy-error-pressure.txt", "w")
    f2 = open("output/lockin_toy-error-diameter.txt", "w")
    f3 = open("output/lockin_toy-error-temperature.txt", "w")
    print("# pressure[bar]","error[%]",file=f1)
    print("# diameter[m]","error[%]",file=f2)
    print("# temperature[K]","error[%]",file=f3)

    # Parameters used
    print("\nEnergy     : ",energy, " eV")
    #print("Pressure:    ",pressure, "bar")
    #print("T/Tc:        ",temperature/temperature_critical_superfluid(pressure))

    global error
    global value
    error = np.array([])
    value = np.array([])

    # Starts the toy simulation for a range of PRESSURES, fixed T/Tc and diameter
    print("\nPRESSURE...")
    diameter = 200e-9;      # [m] vibrating wire diameter
    ttc=0.1  # T/Tc

    Run_Toy_Pressure(1, 30, 1, unused, ttc, diameter, f1)

    # Plot results
    plt.title(str(diameter*1e9)+' nm - T/Tc='+str(ttc))
    plt.plot(value,error, linestyle='', marker='o', color="black")
    plt.xlabel('Pressure [bar]')
    plt.ylabel('Error [%]')
    plt.yscale('log')
    plt.savefig('lockin-error-pressure.pdf')
    plt.savefig('lockin-error-pressure.png')
    plt.show()

    f1.close()

    
    # Starts the toy simulation for a range of DIAMETERS, fixed T/Tc and pressure
    print("\nDIAMETERS...")
    pressure = 5     # [bar] pressure
    ttc=0.1  # T/Tc

    error = np.array([])
    value = np.array([])
    Run_Toy_Diameter(50e-9, 1000e-9, 100e-9, pressure, ttc, unused, f2)

    # Plot results
    plt.title(str(pressure)+' bar - T/Tc='+str(ttc))
    plt.plot(value*1e9,error, linestyle='', marker='o', color="black")
    plt.xlabel('diameter [nm]')
    plt.ylabel('Error [%]')
    plt.yscale('log')
    plt.savefig('lockin-error-diameter.pdf')
    plt.savefig('lockin-error-diameter.png')
    plt.show()

    f2.close()    


    # Starts the toy simulation for a range of T/Tc
    print("\nT/Tc...")
    pressure = 5     # [bar] pressure
    diameter = 200e-9;      # [m] vibrating wire diameter
    
    error = np.array([])
    value = np.array([])
    Run_Toy_Temperature(0.1, 0.18, 0.01, pressure, unused, diameter, f3)

    # Plot results
    plt.title(str(pressure)+' bar - '+str(diameter*1e9)+' nm')
    plt.plot(value,error, linestyle='', marker='o', color="black")
    plt.xlabel('T/T$_c$')
    plt.ylabel('Error [%]')
    plt.yscale('log')
    plt.savefig('lockin-error-temperature.pdf')
    plt.savefig('lockin-error-temperature.png')
    plt.show()

    f3.close()
    
    

