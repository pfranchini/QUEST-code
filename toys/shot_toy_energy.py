'''
Toy MC simulation for the bolometer response
using the QP SHOT noise

Input:
 - Energy

 - Helium-3 Pressure: pressure
 - Base temperature: t_base
 - Wire diameter: d
 - Wire length: l
 - Reponse time; t_w
 - Decay constant: t_b

Output:
 - Error on Energy 

P. Franchini, 1/2023

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.stats import norm

# import Tsepelin code
exec(open("../mod_helium3.py").read())

# import QUEST-DMC code
exec(open("../mod_quest.py").read())

# Configuration file with all the parameters                                                                                                                               
exec(open("config.py").read())

# Plotting style
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'lines.linewidth': 3})
plt.rcParams.update({'xtick.direction': 'in'})
plt.rcParams.update({'ytick.direction': 'in'})

## Parameters ################################################

#volume = 1e-6      # [m^3] Helium-3 cell
#density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)   

#=============================================================

#pressure = 0     # [bar] pressure
#ttc=0.10
#t_base = 150e-6 # [K] base temperature
#d = 400e-9;      # [m] vibrating wire diameter
#l = 2e-3         # [m] vibrating wire lenght

#=============================================================

#t_b = 5.00  # [s] decay constant
#t_w = 0.77  # [s] response time - IS NOT CONSTANT -

#=============================================================

N = 100  # number of toys
verbose=False #False # verbosity for plotting

###########################################################


# Define the noise function for shot-noise
def noise(_deltaf,_fb,_pressure,_temperature):
    #bandwidth = np.pi*_fb/2 # Samuli docet
    bandwidth = min(np.pi*_fb/2, lockin_bandwidth) # Samuli docet
    gap = energy_gap_in_low_temp_limit(_pressure)
    mass=mass_effective(_pressure)*atomic_mass # effective mass [kg]
    particle_density=1/(np.pi**2)*np.power(2*mass/Plankbar_const**2,3/2)*np.sqrt(Fermi_momentum(_pressure)**2/(2*mass))*Boltzmann_const*_temperature*np.exp(-gap/(Boltzmann_const*_temperature)) # QP number density, eq.31
    sigma_n=np.sqrt((particle_density*Fermi_velocity(_pressure)*l*diameter/2)/2) # shot-noise, eq.32
    noise = _deltaf/sigma_n * np.sqrt(bandwidth)    # relative error due to the shot-noise in a bandwith???, eq.35
    #print(_deltaf,1/sigma_n,noise)
    return noise

# Define signal fit function Dw vs time
def df(_t,_fb,_d): # time, base width, delta (delta width)
    _t_w = 1/(np.pi*_fb)
    _t1=5
    return _fb + np.heaviside(_t-_t1,1)*(_d*np.power(t_b/_t_w,_t_w/(t_b-_t_w))*(t_b/(t_b-_t_w))*(np.exp(-(_t-_t1)/t_b) - np.exp(-(_t-_t1)/_t_w)))

###########################################################

def Toy(energy):

    print()
    print("Energy:      ",str(energy), " eV")
    # Input delta(width) from input energy
    delta, _ = DeltaWidth_from_Energy(energy,pressure,t_base,diameter)
    
    # Base width from the input base temperature
    f_base = Width_from_Temperature(t_base,pressure,diameter)

    # Response time
    t_w = 1/(np.pi*f_base)
    
    print("Base width:      ",f_base*1000, " mHz")
    print("Width variation: ",delta*1000,  " mHz")
    print("t_w: ",t_w, "s") 
          
    #t = np.linspace(0, 50, 200) # time
    t = np.linspace(4.5, 50, 1000) # time

    base_toy = np.array([])  # base width distribution
    delta_toy = np.array([]) # delta width distribution

    # Repeat the fit N times
    for i in range(N):

        # Delta f vs time
        deltaf = f_base + np.heaviside(t-5,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-5)/t_b) - np.exp(-(t-5)/t_w)))

        # Add noise based on QP shot noise
        for j in range(len(deltaf)):
            deltaf[j] = np.random.normal(deltaf[j],noise(deltaf[j],f_base,pressure,t_base), 1)
            
        # Random noise    
        #noise = np.random.normal(0, 1e-3, len(t))
        #deltaf = deltaf + noise
        #plt.plot(t, deltaf)
        #plt.show()
        
        # Fit the noise'd distribution        
        popt, pcov = curve_fit(df,t,deltaf)
        # errors on base width and width increase
        base_fit, delta_fit = popt

        delta_toy = np.append(delta_toy,delta_fit)
        base_toy = np.append(base_toy,base_fit)

        
    if verbose and energy==10:
        # Plot deltaf(t)
        plt.plot(t, deltaf*1e3,linestyle='',marker='.', color="black")
        plt.plot(t, df(t,*popt)*1e3)
        plt.xlabel('time [s]')
        plt.ylabel('$\Delta f$ [mHz]')
        plt.savefig('deltaf_toy-example'+str(diameter*1e9)+'.pdf')
        plt.show()

    if verbose and energy==1000:
        # Plot voltage(t)
        plt.plot(t, v_h*f_base/deltaf*1e9,linestyle='',marker='.', color="black")
        plt.plot(t, v_h*f_base/df(t,*popt)*1e9)
        plt.xlabel('time [s]')
        plt.ylabel('$V_H$ [nV]')
        plt.savefig('voltage_toy-example'+str(diameter*1e9)+'.pdf')
        plt.show()

    # Plot toy energy distribution
    if verbose and energy==10000:
        plt.hist(delta_toy*alpha/1000,100)
        plt.xlabel('Energy [KeV]')
        plt.ylabel('Entries')
        plt.savefig('energy_distribution'+str(diameter*1e9)+'.pdf')
        plt.show()

    # Plot toy base distribution
    if verbose and energy==10000:
        plt.hist(base_toy,100)
        plt.xlabel('Base width [Hz]')
        plt.ylabel('Entries')
        plt.savefig('base_width_distribution'+str(diameter*1e9)+'.pdf')
        plt.show()


    ## Gaussian fit for base and delta toy distributions
    (delta_mu, delta_sigma) = norm.fit(delta_toy)
    (base_mu, base_sigma) = norm.fit(base_toy)

    print("Fitted base: ",base_mu," ",base_sigma)
    print("Fitted delta: ",delta_mu," ",delta_sigma)

    print(alpha_prime*delta_mu*base_sigma)
    print(alpha*delta_sigma)
    
    energy_error=  np.sqrt( np.power(alpha_prime*delta_mu*base_sigma,2) + np.power(alpha*delta_sigma,2) )
#    energy_error=  np.sqrt( np.power(alpha_prime*delta_mu*base_sigma,2) )

    return energy_error


def Run_Toy(start_energy, end_energy, step):

    global error
    global e

    for energy in np.arange(start_energy, end_energy, step):

        sigma_energy = Toy(energy)
        print(energy, sigma_energy, sigma_energy/energy*100,"%")

        error = np.append(error,sigma_energy/energy*100)
        e = np.append(e,energy)

        print(energy,*(sigma_energy/energy*100),file=f)


###########################################################


if __name__ == "__main__":

   
    #d=4.5e-6
    _pressure=0
    _temperature=ttc*temperature_critical_superfluid(pressure)
    gap = energy_gap_in_low_temp_limit(_pressure)
    mass=mass_effective(_pressure)*atomic_mass
    flux=1/(np.pi**2)*np.power(2*mass/Plankbar_const**2,3/2)*np.sqrt(Fermi_momentum(_pressure)**2/(2*mass))*Boltzmann_const*_temperature*np.exp(-gap/(Boltzmann_const*_temperature))

    sigma_n=np.sqrt((flux*Fermi_velocity(_pressure)*l*diameter/2)/2)

    Vm=36.8*1e-6   #molar volume, cm^ 3/mol, from Slazav Library
    n1=6e23/Vm    #number density in SI units calculated with avogadro number
    
    print("Effective mass coeff.: ", mass_effective(_pressure))
    print("Flux: ",flux)
    print("He3 density",n1)
    print("Flux/He_density: ",flux/n1)
    
    print("fractional noise: [%]",1/sigma_n*100)

    sigma = np.array([])
    temperature = np.array([])
    for _temperature in np.arange(70e-6, 200e-6, 10e-6):
        particle_density=1/(np.pi**2)*np.power(2*mass/Plankbar_const**2,3/2)*np.sqrt(Fermi_momentum(_pressure)**2/(2*mass))*Boltzmann_const*_temperature*np.exp(-gap/(Boltzmann_const*_temperature)) # QP number density, eq.31
        sigma_n=np.sqrt( (particle_density*Fermi_velocity(_pressure)*l*diameter/2)/2 ) # shot-noise, eq.32
        sigma=np.append(sigma,1/sigma_n) # fractional noise
        temperature=np.append(temperature,_temperature)

    plt.xlabel('Temperature [uK]')
    plt.ylabel('fractional noise')
    plt.plot(temperature*1e6, sigma)
    plt.yscale('log')
    plt.show()


    # Base temperature defined from the T/Tc
    t_base=ttc*temperature_critical_superfluid(pressure)

    # Output file
    f = open("output/shot-error.txt", "w")
    print("# energy[ev]","error[%]",file=f)

    # Parameters used
    print()
    print("Temperature: ",t_base*1e6, " uk")
    print("Diameter:    ",diameter*1e9," nm")
    print("Pressure:    ",pressure, "bar")

    # Calculates alpha_prime for the error propagation
    global alpha_prime
    epsilon=1e-9
    _, alpha1 = DeltaWidth_from_Energy(1000, pressure, t_base - epsilon/2, diameter)
    _, alpha2 = DeltaWidth_from_Energy(1000, pressure, t_base + epsilon/2, diameter)
    gap = energy_gap_in_low_temp_limit(pressure)
    alpha_prime = (alpha2-alpha1)/epsilon * Boltzmann_const/gap * np.power(t_base,2) * 1/Width_from_Temperature(t_base,pressure,diameter)
    print("alpha_prime",alpha_prime)

    # Single toy run for a defined energy [eV]
    #_ = Toy(10);
    
    # Starts the toy simulation for a range of energies
    global error
    global e
    error = np.array([])
    e = np.array([])

    Run_Toy(1e-2, 0.9, 1e-1)
    Run_Toy(1, 9, 1)
    Run_Toy(10, 90, 10)
    Run_Toy(100, 900, 50)
    Run_Toy(1000, 9000, 500)
    Run_Toy(10000, 100000, 2000)

    # Close file
    f.close()
    
    # Plot results
    plt.title(str(diameter*1e9)+" nm - "+str(t_base*1e6)+" $\mu$K - "+str(pressure)+ " bar")
    plt.plot(e/1000,error, linestyle='', marker='o', color="black")
    plt.xscale('log')
    plt.ylim([0, 100])  
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Error [%]')
    plt.savefig('error-shot'+str(diameter*1e9)+'.png')
    plt.show()


