'''
Toy MC simulation for the bolometer response
using the error from the SQUID circuit

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
exec(open("mod_helium3.py").read())

## Parameters ################################################

volume = 1e-6      # [m^3] Helium-3 cell
density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)   

#=============================================================

pressure = 5     # [bar] pressure
t_base = 150e-6  # [K] base temperature
d = 200e-9;      # [m] vibrating wire diameter

#=============================================================

t_b = 5.00  # [s] decay constant
#t_w = 0.77  # [s] response time - IS NOT CONSTANT -

v_h = np.pi/2*1e-7  # [V] Base voltage height for a v=1mm/s

#=============================================================

N = 100  # number of toys
verbose=False # verbosity for plotting

## SQUID parameters ==========================================                                                                                                             

B = 0.4e-3 # T
R = 1 # Ohm
w0 = 5000  # Hz
L = 1.5e-6 # H
#Z = w0*L  minimum at fields
phi0 = 2.067833848e-15  # Wb
S = np.power(0.4e-6*phi0,2) # Hz^-1
M = 10.0e-9 # H
l = 2e-3 # m
m = np.power(d/2,2)*np.pi*density # kg/m
v0 = 1e-3 # m/sec

# ===========================================================

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


# Define signal fit function Dw vs time
def df(_t,_fb,_d): # time, base width, delta (delta width)
    _t_w = 1/(np.pi*_fb)
    _t1=5
    return _fb + np.heaviside(_t-_t1,1)*(_d*np.power(t_b/_t_w,_t_w/(t_b-_t_w))*(t_b/(t_b-_t_w))*(np.exp(-(_t-_t1)/t_b) - np.exp(-(_t-_t1)/_t_w)))

def Voltage_Error(_fb): # base width
    Z=l*np.power(B,2)/(2*np.pi*m*_fb)
    v_rms = np.sqrt( (np.power(Z+R,2) + np.power(w0*L,2))*S*np.pi*_fb/2/np.power(M,2) + 4*Boltzmann_const*t_base*R*np.pi*_fb/2 + Boltzmann_const*t_base*l*np.power(B,2)/m ) # [V]
    nep = 0.250 # Noise Equivalent Power in units of bandwidth
    return v_rms*np.sqrt(nep)


###########################################################

def Toy(energy):

    print()
    print("Energy:      ",str(energy), " eV")
    # Input delta(width) from input energy
    delta, _ = DeltaWidth_from_Energy(energy,pressure,t_base)
    
    # Base width from the input base temperature
    f_base = Width_from_Temperature(t_base,pressure)

    # Response time
    t_w = 1/(np.pi*f_base)
    
    print("Base width:      ",f_base*1000, " mHz")
    print("Width variation: ",delta*1000,  " mHz")
    print("t_w: ",t_w, "s")
    print("Bandwidth: ",np.pi*f_base/2, " Hz")
    print("Voltage error:",Voltage_Error(f_base))

    
    t = np.linspace(0, 50, 200) # time

    base_toy = np.array([])  # base width distribution
    delta_toy = np.array([]) # delta width distribution

    # Repeat the fit N times
    for i in range(N):

        # Delta f vs time
        deltaf = f_base + np.heaviside(t-5,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-5)/t_b) - np.exp(-(t-5)/t_w)))

        # Add noise based on voltage error. The voltage error comes from the SQUID circuit.
        for j in range(len(deltaf)):
            deltaf[j] = np.random.normal(deltaf[j],np.power(deltaf[j],2)/(v_h*f_base)*Voltage_Error(f_base), 1)

        # Fit the noise'd distribution        
        popt, pcov = curve_fit(df,t,deltaf)
        # errors on base width and width increase
        base_fit, delta_fit = popt

        delta_toy = np.append(delta_toy,delta_fit)
        base_toy = np.append(base_toy,base_fit)

        
    if verbose and energy==10000:
        # Plot deltaf(t)
        plt.plot(t, deltaf,linestyle='',marker='.', color="black")
        plt.plot(t, df(t,*popt))
        plt.xlabel('time [s]')
        plt.ylabel('$\Delta f$ [Hz]')
        plt.savefig('deltaf_toy-example'+str(d*1e9)+'.pdf')
        plt.show()

    if verbose and energy==10000:
        # Plot voltage(t)
        plt.plot(t, v_h*f_base/deltaf*1e9,linestyle='',marker='.', color="black")
        plt.plot(t, v_h*f_base/df(t,*popt)*1e9)
        plt.xlabel('time [s]')
        plt.ylabel('$V_H$ [nV]')
        plt.savefig('voltage_toy-example'+str(d*1e9)+'.pdf')
        plt.show()

    # Plot toy energy distribution
    if verbose and energy==10000:
        plt.hist(delta_toy*alpha/1000,100)
        plt.xlabel('Energy [KeV]')
        plt.ylabel('Entries')
        plt.savefig('energy_distribution'+str(d*1e9)+'.pdf')
        plt.show()

    # Plot toy base distribution
    if verbose and energy==10000:
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
        print("Error:", sigma_energy,"eV, ", sigma_energy/energy*100,"%")

        error = np.append(error,sigma_energy/energy*100)
        e = np.append(e,energy)

        print(energy,*(sigma_energy/energy*100),file=f)


###########################################################


if __name__ == "__main__":

    # Output file
    f = open("squid_toy-error.txt", "w")
    print("# energy[ev]","error[%]",file=f)

    # Parameters used
    print()
    print("Temperature: ",t_base*1e6, " uk")
    print("Diameter:    ",d*1e9," nm")
    print("Pressure:    ",pressure, "bar")
    print("T/Tc:        ",t_base/temperature_critical_superfluid(pressure))

    # Calculates alpha_prime for the error propagation
    global alpha_prime
    epsilon=1e-9
    _, alpha1 = DeltaWidth_from_Energy(1000, pressure, t_base - epsilon/2)
    _, alpha2 = DeltaWidth_from_Energy(1000, pressure, t_base + epsilon/2)
    gap = energy_gap_in_low_temp_limit(pressure)
    alpha_prime = (alpha2-alpha1)/epsilon * Boltzmann_const/gap * np.power(t_base,2) * 1/Width_from_Temperature(t_base,pressure)
    print("alpha_prime",alpha_prime)
    
    # Starts the toy simulation for a range of energies
    global error
    global e
    error = np.array([])
    e = np.array([])

    Run_Toy(1e-2, 0.9, 1e-2)
    Run_Toy(1, 9, 1)
    Run_Toy(10, 90, 10) 
    Run_Toy(100, 900, 50)
    Run_Toy(1000, 9000, 500)
    Run_Toy(10000, 100000, 2000)

    f.close()
    
    # Plot results
    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - "+str(t_base*1e6)+" $\mu$K - "+str(pressure)+ " bar")
    plt.plot(e/1000,error, linestyle='', marker='o', color="black")
    plt.xscale('log')
    #plt.ylim([0, 3])  
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Error [%]')
    plt.savefig('error'+str(d*1e9)+'.pdf')
    plt.savefig('error'+str(d*1e9)+'.png')
    plt.show()


    # Single toy run for a defined energy
    #_ = Toy(10000);

    

