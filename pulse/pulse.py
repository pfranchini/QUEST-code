'''

Study of the position of the maximum of the pulse position in the time window.
The time corresponds to peak time, the maximum in the energy deposition pulse,
according to Winkelmann equation pulse shape.

Usage:
    pulse.py -p <pressure [bar]> -d <diameter [m]> -t <relative temperature> -b <decay constant t_b>

'''

import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# import Tsepelin code
exec(open("../mod_helium3.py").read())

#=============================================================  

density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)

#=============================================================

def arguments(argv):

    global arg_pressure
    global arg_diameter
    global arg_ttc
    global arg_t_b
    
    arg_pressure = 5
    arg_diameter = 200e-9
    arg_ttc = 0.1
    arg_t_b = 5.00 # [s] decay constant  
    
    arg_help = "{0} -p <pressure [bar]> -d <diameter [m]> -t <relative temperature> -b <decay constant t_b>".format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hp:d:t:b:", ["help", "pressure=", "diameter=", "ttc=", "t_b="])
    except:
        print(arg_help)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            
            sys.exit(2)
        elif opt in ("-p", "--pressure"):
            arg_pressure = float(arg)
        elif opt in ("-d", "--diameter"):
            arg_diameter = float(arg)
        elif opt in ("-t", "--ttc"):
            arg_ttc = float(arg)
        elif opt in ("-b", "--t_b"):
            arg_t_b = float(arg)

#=============================================================
            
def Width_from_Temperature(_temperature,_pressure,_diameter):
    gap = energy_gap_in_low_temp_limit(_pressure)
    width=np.power(Fermi_momentum(_pressure),2)*Fermi_velocity(_pressure)*density_of_states(_pressure)/(2*density*np.pi*_diameter)*np.exp(-gap/(Boltzmann_const*_temperature))
    return width


def t_prime(_fb):
    t_w = 1/(np.pi*_fb)
    return t_b*t_w/(t_w-t_b)*np.log(t_w/t_b)

#=============================================================

if __name__ == "__main__":

    arguments(sys.argv)
    t_b=arg_t_b
    
    # Parameters used
    tmin=0.1 # [s]
    tmax=10  # [s]
    
    print("Parameters:")
    print("Diameter:    ",arg_diameter*1e9," nm")
    print("Pressure:    ",arg_pressure, " bar")
    print("T/Tc:        ",arg_ttc)
    print("t_b:         ",arg_t_b, " s")

    
    # t_prime vs DIAMETER
    pressure = arg_pressure
    ttc=arg_ttc

    diameter = np.array([])
    time = np.array([])

    for d in np.arange(50e-9, 1000e-9, 100e-9):
        time =  np.append(time, t_prime(Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,d)))
        diameter = np.append(diameter, d)

    # Plot results
    plt.title(str(pressure)+' bar - T/Tc='+str(ttc))
    plt.plot(diameter*1e9,time, linestyle='', marker='o', color="black")
    plt.axhline(tmin,color='grey', linestyle='--')
    plt.axhline(tmax,color='grey', linestyle='--')
    plt.xlabel('diameter [nm]')
    plt.ylabel('peak time [s]')
    plt.savefig('pulse_diameter.png')
    plt.show()


    # t_prime vs PRESSURE
    ttc=arg_ttc
    diameter=arg_diameter
    
    pressure = np.array([])
    time = np.array([])

    for p in np.arange(1, 30, 1,):
        time =  np.append(time, t_prime(Width_from_Temperature(ttc*temperature_critical_superfluid(p),p,diameter)))
        pressure = np.append(pressure, p)

    # Plot results
    plt.title(str(diameter*1e9)+' nm - T/Tc='+str(ttc))
    plt.plot(pressure,time, linestyle='', marker='o', color="black")
    plt.axhline(tmin,color='grey', linestyle='--')
    plt.axhline(tmax,color='grey', linestyle='--')
    plt.xlabel('pressure [bar]')
    plt.ylabel('peak time [s]')
    plt.savefig('pulse_pressure.png')
    plt.show()


    # t_prime vs TEMPERATURE
    pressure = arg_pressure
    diameter = arg_diameter
    
    TTc = np.array([])
    time = np.array([])

    for ttc in np.arange(0.1, 0.3, 0.01):
        time = np.append(time, t_prime(Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,diameter)))
        TTc = np.append(TTc, ttc)

    # Plot results
    plt.title(str(pressure)+' bar - '+str(diameter*1e9)+' nm')
    plt.plot(TTc,time, linestyle='', marker='o', color="black")
    plt.axhline(tmin,color='grey', linestyle='--')
    plt.axhline(tmax,color='grey', linestyle='--')
    plt.xlabel('T/T$_c$')
    plt.ylabel('peak time [s]')
    plt.savefig('pulse_temperature.png')
    plt.show()
