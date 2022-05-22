'''
Analysis of Viktor's data
=========================

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
 - Energy spectra

'''

import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.signal import find_peaks

# import Tsepelin code
exec(open("mod_helium3.py").read())

## Input #####################################################

energy = 10000 # [eV]


## Parameters ################################################

volume = 1.25e-7    # [m^3] Black box radiator
density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)   

#=============================================================

pressure = 0    # [bar] pressure
#t_base = 111e-6  # [K] base temperature
d = 4.5e-6;      # [m] 3m2: vibrating wire diameter
#d = 13.5e-6;      # [m] 3m2: vibrating wire diameter

t_b = 5.00  # [s] decay constant
#t_w = 0.77  # [s] response time - IS NOT CONSTANT -

v_h = np.pi/2*1e-7  # [V] Base voltage height for a v=1mm/s
v_rms = 3.5*1e-9    # [V] Error on voltage measurement for a lock-in amplifier
# v_drive=14.5e-3

#=============================================================

N = 1000  # number of toys
verbose=True # verbosity for plotting

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
    if verbose:
        plt.plot(DQ/1e3,DW*1e3,label='DQvsDW')
        plt.title('Width variation vs Deposited energy')
        plt.xlim([0, 100])
        plt.ylim([0, 200])
        plt.xlabel('$\Delta$Q [KeV]')
        plt.ylabel('$\Delta(\Delta f)$ [mHz]')
        plt.show()
    
    
    # Fit line to extract the slope alpha: DQ = alpha * DW
    global alpha
    alpha, _ = np.polyfit(DW, DQ, 1)
    
    # Input delta width from input energy
    deltawidth = E/alpha
       
    return deltawidth, alpha

# Define signal fit function Dw vs time
#def df(_t,_fb,_d,_t_w,_t1,t_b): # time, base width, delta (delta width)
#    _t_w = 1/(np.pi*_fb)
#    _t1=5
#    return _fb + np.heaviside(_t-_t1,1)*(_d*np.power(t_b/_t_w,_t_w/(t_b-_t_w))*(t_b/(t_b-_t_w))*(np.exp(-(_t-_t1)/t_b) - np.exp(-(_t-_t1)/_t_w)))

def df(time,a,b):
    return a*time+b


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
    
    t = np.linspace(0, 50, 200) # time

    base_toy = np.array([])  # base width distribution
    delta_toy = np.array([]) # delta width distribution

    # Repeat the fit N times
    for i in range(N):

        # Delta f vs time
        deltaf = f_base + np.heaviside(t-5,1)*(delta*np.power(t_b/t_w,t_w/(t_b-t_w))*(t_b/(t_b - t_w))*(np.exp(-(t-5)/t_b) - np.exp(-(t-5)/t_w)))

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
    
###########################################################


if __name__ == "__main__":

    print()
    #print("Temperature: ",t_base*1e6, " uk")
    print("Diameter:    ",d*1e9," nm")
    print("Pressure:    ",pressure, "bar")


    # Reads the data
    #with open('/home/franchini/Documents/QUEST/Fridge5_AERO5_Run4_DM13_TWDAQ05_analysed-REDUCED.dat', 'rt') as f:
    with open('/home/franchini/Documents/QUEST/Fridge5_AERO5_Run4_DM13_TWDAQ05_analysed.dat', 'rt') as f:

        reader = csv.reader(f, delimiter=',', skipinitialspace=True)
        lineData = list()
        cols = next(reader)
  
        # Print names of variables
        print(cols)
    
        for col in cols:
            # Create a list in lineData for each column of data.
            lineData.append(list())
        
        for line in reader:
            for i in range(0, len(lineData)):
                # Copy the data from the line into the correct columns.
                lineData[i].append(line[i])
                
        data = dict()
    
        for i in range(0, len(cols)):
            # Create each key in the dict with the data in its column.
            data[cols[i]] = lineData[i]
            
    #=======================================================================

    t = np.array(data['secs'])
    width_ = np.array(data['TW27width'])
    temperature_ = np.array(data['TW27temperature'])

    width = np.array([])
    for w in width_:
        width = np.append(width,float(w))

    temperature = np.array([])
    for temp in temperature_:
        temperature = np.append(temperature,float(temp))

    print(temperature)
        
    t_base = np.nanmean(temperature)
    print("Average temperature:      ",t_base, " K")

    plt.xlabel('time [s]')
    plt.ylabel('temperature [uK]')
    plt.plot(t,temperature*1e6)
    plt.show()
  
    # Base width from the input base temperature                                                                                                                           
    f_base = Width_from_Temperature(t_base,pressure)
    print("Calculated base width:      ",f_base*1000, " mHz")

    background = np.nanmean(width)
    print("Background: ",background*1000, " mHz")

    plt.xlabel('time [s]')
    plt.ylabel('width [Hz]')
    plt.plot(t,width)
    plt.show()

    threshold = 0 #0.3 #0.2
    peaks, _ = find_peaks(width, height=(f_base[0]+threshold))

    print("NUmber of peaks: ",len(peaks))
    
    _, alpha = DeltaWidth_from_Energy(1000,pressure,t_base)
    print("Alpha: ",alpha)

    energy = (width[peaks]-background)*alpha

    plt.plot(width)
    plt.plot(peaks, width[peaks], "x")
    plt.plot(np.zeros_like(width)+f_base[0]+threshold, "--", color="gray")
    plt.show()

    plt.hist(energy/1e3, 400, range=[0, 2000])
    plt.yscale('log')
    plt.xlim([0, 2000])    
    plt.xlabel('Energy [keV]')
    plt.show()     

    '''
    #interval = (x>=9821) & (x<=9830)
    plt.xlim([9821,9830])
    plt.xlabel('time [s]')
    plt.ylabel('width [Hz]')
    plt.plot(t,width)
    plt.show()

    #width2 = np.array([])
    #for w in width:
    #    width2 = np.append(width2,float(w))

    #plt.hist(width2, 100)
    #plt.xlim([0,0.8])
    #plt.show()
    #    plt.plot(x,y[(x>=9821) & (x<=9830)])
    #plt.plot(data['secs'],data['TW29temperature'])
    '''




#    popt, pcov = curve_fit(df,x,y)

    '''
    # Calculates alpha_prime for the error propagation
    global alpha_prime
    epsilon=1e-9
    _, alpha1 = DeltaWidth_from_Energy(energy, pressure, t_base - epsilon/2)
    _, alpha2 = DeltaWidth_from_Energy(energy, pressure, t_base + epsilon/2)
    gap = energy_gap_in_low_temp_limit(pressure)
    alpha_prime = (alpha2-alpha1)/epsilon * Boltzmann_const/gap * np.power(t_base,2) * 1/Width_from_Temperature(t_base,pressure)

    print("alpha_prime",alpha_prime)

    
    # Starts the toy simulation for a range of energies
    error = np.array([])
    e = np.array([])

    
    for energy in np.arange(100, 900, 50):

        sigma_energy = Toy(energy)
        print(energy, sigma_energy, sigma_energy/energy*100,"%")

        error = np.append(error,sigma_energy/energy*100)
        e = np.append(e,energy)

    for energy in np.arange(1000, 9000, 500):

        sigma_energy = Toy(energy)
        print(energy, sigma_energy, sigma_energy/energy*100,"%")

        error = np.append(error,sigma_energy/energy*100)
        e = np.append(e,energy)


    for energy in np.arange(10000, 100000, 2000):

        sigma_energy = Toy(energy)
        print(energy, sigma_energy, sigma_energy/energy*100,"%")

        error = np.append(error,sigma_energy/energy*100)
        e = np.append(e,energy)


    
    plt.plot(e/1000,error, linestyle='', marker='o', color="black")
    plt.xscale('log')
    plt.ylim([0, 100])  
    plt.xlabel('Energy [KeV]')
    plt.ylabel('Error [%]')
    plt.savefig('error'+str(d*1e9)+'.pdf')
    plt.savefig('error'+str(d*1e9)+'.png')
    plt.show()

    
            

    _ = Toy(10000);

    '''    
