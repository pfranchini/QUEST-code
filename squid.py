'''

Simulation for the bolometer response with a SQUID readout
considering the resolution as in Lev's note.

-----------------------------------------------------------

Usage:
  squid.py -p <pressure [bar]> -d <diameter [m]> -t <relative temperature> -B <magnetic field>

-----------------------------------------------------------

Input:
 - Energy

 - Helium-3 Pressure: pressure
 - Wire diameter: d
 - Base temperature: t_base
 - central resonating frequency

 - SQUID parameters

Output:
 - Resolution vs Pressure and T/Tc for each noise component
 - Resolution vs Pressure, with R dependency
 - Resolution vs B, with R dependency
 - Error on Energy, with R dependency


'''

import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.stats import norm

# import Tsepelin code
exec(open("mod_helium3.py").read())


## Parameters ################################################

volume = 1e-6      # [m^3] Helium-3 cell
density = 6.05e3  # [kg/m^3] Niobium-Titanium (NbTi)   

#=============================================================

pressure = 5     # [bar] pressure
TTc = 0.1        # relative temperature
d = 200e-9       # [m] vibrating wire diameter
B = 0.4e-3       # [T] magnetic field

#=============================================================

verbose=False # verbosity for plotting

## SQUID parameters ==========================================

#B = 0.4e-3 # T # as field uses always the one that minimize dW
R = 1 # Ohm

w0 = 5000  # Hz
L = 1.5e-6 # H
#Z = w0*L  minimum at fields
phi0 = 2.067833848e-15  # Wb
S = np.power(0.4e-6*phi0,2) # Hz^-1
M = 10.0e-9 # H
l = 2e-3 # m
m = np.power(d/2,2)*np.pi*density # kg/m

# ===========================================================
   

## More routines: ###########################################

def arguments(argv):

    global pressure
    global d
    global TTc
    global B
    arg_help = "{0} -p <pressure [bar]> -d <diameter [m]> -t <relative temperature> -B <magnetic field>".format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hp:d:t:B:", ["help", "pressure=", "diameter=", "ttc=", "field="])
    except:
        print(arg_help)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-p", "--pressure"):
            pressure = float(arg)
        elif opt in ("-d", "--diameter"):
            d = float(arg)
        elif opt in ("-t", "--ttc"):
            TTc = float(arg)
        elif opt in ("-B", "--field"):
            B = float(arg)
            

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
def SQUID_Resolution(_pressure,_temperature,magnetic_field,resistance):
    f_base = Width_from_Temperature(_temperature,_pressure)
    Z=l*np.power(magnetic_field,2)/(2*np.pi*m*f_base)
    v0=Boltzmann_const*_temperature/Fermi_momentum(_pressure)
    V=l*v0*magnetic_field # max voltage across the wire

    return np.sqrt( (np.power(Z+resistance,2) + np.power(w0*L,2))*S*np.pi*f_base/2/np.power(M,2) + 4*Boltzmann_const*_temperature*resistance*np.pi*f_base/2 + Boltzmann_const*_temperature*l*np.power(magnetic_field,2)/m )/V * f_base


# SQUID resolution component as in Lev's note (dW) [Hz]
def SQUID_Resolution_A(_pressure,_temperature,magnetic_field,resistance):
    f_base = Width_from_Temperature(_temperature,_pressure)
    Z=l*np.power(magnetic_field,2)/(2*np.pi*m*f_base)
    v0=Boltzmann_const*_temperature/Fermi_momentum(_pressure)
    V=l*v0*magnetic_field # max voltage across the wire

    return np.sqrt( (np.power(Z+resistance,2) + np.power(w0*L,2))*S*np.pi*f_base/2/np.power(M,2) )/V * f_base


# SQUID resolution component as in Lev's note (dW) [Hz]
def SQUID_Resolution_B(_pressure,_temperature,magnetic_field,resistance):
    f_base = Width_from_Temperature(_temperature,_pressure)
    Z=l*np.power(magnetic_field,2)/(2*np.pi*m*f_base)
    v0=Boltzmann_const*_temperature/Fermi_momentum(_pressure)
    V=l*v0*magnetic_field # max voltage across the wire

    return np.sqrt( 4*Boltzmann_const*_temperature*resistance*np.pi*f_base/2 )/V * f_base


# SQUID resolution component as in Lev's note (dW) [Hz]
def SQUID_Resolution_C(_pressure,_temperature,magnetic_field,resistance):
    f_base = Width_from_Temperature(_temperature,_pressure)
    Z=l*np.power(magnetic_field,2)/(2*np.pi*m*f_base)
    v0=Boltzmann_const*_temperature/Fermi_momentum(_pressure)
    V=l*v0*magnetic_field # max voltage across the wire

    return np.sqrt( Boltzmann_const*_temperature*l*np.power(magnetic_field,2)/m )/V * f_base

###########################################################

# error calculation (is not a toy anymore)
def Toy(energy):
    
    print()
    print("Energy:      ",str(energy), " eV")
    # Input delta(width) from input energy
    delta, _ = DeltaWidth_from_Energy(energy,pressure,TTc*temperature_critical_superfluid(pressure))
    
    # Base width from the input base temperature
    f_base = Width_from_Temperature(TTc*temperature_critical_superfluid(pressure),pressure)

    print("Base width:      ",f_base*1000, " mHz")
    print("Width variation: ",delta*1000,  " mHz")

    # Calculate the resolution from eq.8 in the Lev's note on SQUID
    #for field in np.arange(1e-6,1,1e3): 
    delta_sigma = SQUID_Resolution(pressure,temperature,B,R)
    
    energy_error=  alpha*delta_sigma


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
    f = open("squid_error.txt", "w")
    print("# energy[ev]","error[%]",file=f)

    arguments(sys.argv)

    # Parameters used
    print("Parameters:")
    print("Diameter:    ",d*1e9," nm")
    print("Length:      ",l*1e3, "mm")
    print("Pressure:    ",pressure, "bar")
    print("T/Tc:        ",TTc)
    print("Field:       ",B, "T")
    print("Temperature: ",TTc*temperature_critical_superfluid(pressure)*1e6, " uk")

    energy= 10000
    delta, _ = DeltaWidth_from_Energy(energy,pressure,TTc*temperature_critical_superfluid(pressure))
    
    print("Thermal motion energy: ",alpha*Fermi_momentum(pressure)/np.sqrt(l*m*Boltzmann_const*TTc*temperature_critical_superfluid(pressure))*Width_from_Temperature(TTc*temperature_critical_superfluid(pressure),pressure)
, "eV")

    '''
    # Preliminary validation plots;

    # Base width vs pressure
    Pressure = np.array([])  # [bar]
    Width = np.array([])  # [Hz]
    for p in np.arange(0, 30, 0.1):
        Pressure=np.append(Pressure,p)
        Width=np.append(Width,Width_from_Temperature(t_base,p))

    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - "+str(t_base*1e6)+" $\mu$K")
    plt.plot(Pressure, Width)
    plt.yscale('log')
    plt.xlabel('Pressure [bar]')
    plt.ylabel('Base width [Hz]')
    plt.show()

    # Base width vs pressure for T/Tc constant
    for TTc in (0.10, 0.15, 0.20):
        Pressure = np.array([])  # [bar]
        Width = np.array([])  # [Hz]
        for p in np.arange(0, 30, 0.1):
            temperature=TTc*temperature_critical_superfluid(p)*0.92
            Pressure=np.append(Pressure,p)
            Width=np.append(Width,Width_from_Temperature(temperature,p))

        plt.plot(Pressure, Width,  label='$T/T_c$='+str(TTc))
            
    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm")
    plt.yscale('log')
    plt.xlabel('Pressure [bar]')
    plt.ylabel('Base width [Hz]')
    plt.legend()
    plt.show()

    # Alpha width vs pressure
    Pressure = np.array([])  # [bar]
    Alpha = np.array([])  # [Hz]
    for p in np.arange(0, 30, 1):
        Pressure=np.append(Pressure,p)
        _, a = DeltaWidth_from_Energy(10000,p,t_base)
        Alpha=np.append(Alpha,a)

    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - "+str(t_base*1e6)+" $\mu$K")
    plt.plot(Pressure, Alpha)
    plt.xlabel('Pressure [bar]')
    plt.ylabel('Alpha [eV/Hz]')
    plt.show()
    '''

    # Resolution vs pressure, divided by contributions
    print("Resolution vs pressure...")

    R=0.001
    
    Pressure = np.array([])  # [bar]
    Resolution = np.array([])  # [Hz]
    ResolutionA = np.array([])  # [Hz]
    ResolutionB = np.array([])  # [Hz]
    ResolutionC = np.array([])  # [Hz]

    for p in np.arange(0, 30, 0.1):
        temperature=TTc*temperature_critical_superfluid(p)
        Pressure=np.append(Pressure,p)
        Resolution=np.append(Resolution,SQUID_Resolution(p,temperature,B,R))
        ResolutionA=np.append(ResolutionA,SQUID_Resolution_A(p,temperature,B,R))
        ResolutionB=np.append(ResolutionB,SQUID_Resolution_B(p,temperature,B,R))
        ResolutionC=np.append(ResolutionC,SQUID_Resolution_C(p,temperature,B,R))

    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - $T/T_c$="+str(TTc))
    plt.plot(Pressure, Resolution, label="TOTAL")
    plt.plot(Pressure, ResolutionA, label="current noise")
    plt.plot(Pressure, ResolutionB, label="Johnson-Nyquist noise")
    plt.plot(Pressure, ResolutionC, label="thermal noise")
    #plt.yscale('log')
    plt.xlabel('Pressure [bar]')
    plt.ylabel('Resolution [Hz]')
    #plt.axhline(y=thermal_motion, color='grey', linestyle='-')
    plt.savefig("output/squid/sensitivity_vs_pressure.png")
    plt.legend()
    plt.savefig("output/squid/sensitivity_vs_pressure.png")
    plt.show()


    # Resolution vs T/Tc, divided by contributions
    print("Resolution vs relative temperature...")
    R=0.001
    
    Temperature = np.array([])  # 
    Resolution = np.array([])   # [Hz]
    ResolutionA = np.array([])  # [Hz]
    ResolutionB = np.array([])  # [Hz]
    ResolutionC = np.array([])  # [Hz]

    for ttc in np.arange(0.1, 0.15, 0.01):
        t=ttc*temperature_critical_superfluid(pressure)

        Temperature=np.append(Temperature,ttc)
        Resolution=np.append(Resolution,SQUID_Resolution(pressure,t,B,R))
        ResolutionA=np.append(ResolutionA,SQUID_Resolution_A(pressure,t,B,R))
        ResolutionB=np.append(ResolutionB,SQUID_Resolution_B(pressure,t,B,R))
        ResolutionC=np.append(ResolutionC,SQUID_Resolution_C(pressure,t,B,R))

    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - Pressure="+str(pressure))
    plt.plot(Temperature, Resolution, label="TOTAL")
    plt.plot(Temperature, ResolutionA, label="current noise")
    plt.plot(Temperature, ResolutionB, label="Johnson-Nyquist noise")
    plt.plot(Temperature, ResolutionC, label="thermal noise")
    plt.yscale('log')
    plt.xlabel('T/T$_c$')
    plt.ylabel('Resolution [Hz]')
    #plt.axhline(y=thermal_motion, color='grey', linestyle='-')
    plt.legend()
    plt.savefig("output/squid/sensitivity_vs_temperature.png")
    plt.show()

   
    # Resolution vs pressure, for different wire/contact resistances and T/Tc constant
    print("Resolution vs pressure...")
    for R in (1, 0.1, 0.01, 0.001, 0):

        Pressure = np.array([])  # [bar]
        Resolution = np.array([])  # [Hz]

        for p in np.arange(0, 30, 0.1):
            temperature=TTc*temperature_critical_superfluid(p)
            f_base = Width_from_Temperature(temperature,p)
            Pressure=np.append(Pressure,p)
            Resolution=np.append(Resolution,SQUID_Resolution(p,temperature,B,R))

        plt.plot(Pressure, Resolution, label=str(R)+' $\Omega$')

    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - T/Tc="+str(TTc))
    #plt.yscale('log')
    plt.xlabel('Pressure [bar]')
    plt.ylabel('Resolution [Hz]')
    #plt.axhline(y=thermal_motion, color='grey', linestyle='-')
    plt.legend()
    plt.show()


    # Resolution vs the magnetic field, for different wire/contact resistances
    print("Resolution vs magnetic field...")
    f_base = Width_from_Temperature(TTc*temperature_critical_superfluid(pressure),pressure)
    for R in (1, 0.1, 0.01, 0.001, 0):

        print(" resistance: ",R)

        Field = np.array([])  # [T]
        Resolution = np.array([])  # [Hz]

        for bfield in np.arange(1e-6, 1, 1e-5):
            Field=np.append(Field,bfield)
            Resolution=np.append(Resolution,SQUID_Resolution(pressure,TTc*temperature_critical_superfluid(pressure),bfield,R))

        plt.plot(Field, Resolution, label=str(R)+' $\Omega$')

        # Threshold from dW for field corresponding to the minimum 
        print("  Threshold: ",alpha*SQUID_Resolution(pressure,TTc*temperature_critical_superfluid(p),Field[Resolution.argmin()],R), "eV")

    thermal_motion=Fermi_momentum(pressure)/np.sqrt(l*m*Boltzmann_const*TTc*temperature_critical_superfluid(p))*Width_from_Temperature(TTc*temperature_critical_superfluid(p),pressure)

    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - T/Tc="+str(TTc)+" - "+str(pressure)+ " bar")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Magnetic field [T]')
    plt.ylabel('Resolution [Hz]')
    plt.axhline(y=thermal_motion, color='grey', linestyle='-')
    plt.legend()
    plt.show()


    # Error for different wire/contact resistances
    temperature=TTc*temperature_critical_superfluid(pressure)
    print("Error vs energy...")
    for R in (1, 0.1, 0.01, 0.001, 0):

        print("Resistance: ",R)

        # Finds the field that minimize dW
        Field = np.array([])  # [T]
        Resolution = np.array([])  # [Hz]

        for bfield in np.arange(1e-6, 1, 1e-5):
            Field=np.append(Field,bfield)
            Resolution=np.append(Resolution,SQUID_Resolution(pressure,temperature,bfield,R))

        # Threshold from dW for field corresponding to the minimum 
        B = Field[Resolution.argmin()]

        error = np.array([])
        e = np.array([])

        Run_Toy(1, 100, 10)
        Run_Toy(100, 900, 10)
        Run_Toy(1000, 9000, 500)
        Run_Toy(10000, 100000, 2000)

        # Draw a line for each resistance
        plt.plot(e/1000,error,label=str(R)+' $\Omega$')     

    # Plot results
    plt.title(str(d*1e9)+" nm - "+str(l*1e3)+" mm - $T/T_c$="+str(TTc) + " $\mu$K - "+str(pressure)+ " bar")
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
        
