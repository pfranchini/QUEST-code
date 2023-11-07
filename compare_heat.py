'''
Compare heat capacity for the two equations
 - V&W (eq. 3.81)
 - Winkelmann (eq. 9)
'''

import numpy as np
import matplotlib.pyplot as plt

# import Tsepelin code
exec(open("mod_helium3.py").read())

def heat_capacity_Cv_B_phase_Winkelmann(TemperatureK, PressureBar):

    C0=1.7e3 # [J/K/m³] as in Winkelmann
    gap = energy_gap_in_low_temp_limit(PressureBar)
    cv = C0 * np.power(temperature_critical_superfluid(PressureBar)/TemperatureK, 3./2.) * np.exp(-1 * gap / (Boltzmann_const * TemperatureK))

    return cv

###########################################################

if __name__ == "__main__":

    
    pressure = 0
    T_Tc = np.linspace(0.001, 0.5, 1000)

    print("Winkelmann T/T_c",130e-6/temperature_critical_superfluid(pressure))

    print("T/Tc=0.15:",heat_capacity_Cv_B_phase(0.15*temperature_critical_superfluid(pressure),pressure))
    print("T/Tc=0.17:",heat_capacity_Cv_B_phase(0.17*temperature_critical_superfluid(pressure),pressure))
    print("T/Tc=0.15:",heat_capacity_Cv_B_phase_Winkelmann(0.15*temperature_critical_superfluid(pressure),pressure))

    
    plt.plot(T_Tc, heat_capacity_Cv_B_phase(           T_Tc*temperature_critical_superfluid(pressure),pressure), label='V&W (eq. 3.81)')
    plt.plot(T_Tc, heat_capacity_Cv_B_phase_Winkelmann(T_Tc*temperature_critical_superfluid(pressure),pressure), label='Winkelmann (eq. 9)')
    plt.xlabel('T/$T_c$')
    plt.ylabel('$C_v$ [J/K/m$^3$]')
    plt.legend()
    plt.show()

    
    plt.plot(T_Tc, (heat_capacity_Cv_B_phase(           T_Tc*temperature_critical_superfluid(pressure),pressure)-heat_capacity_Cv_B_phase_Winkelmann(T_Tc*temperature_critical_superfluid(pressure),pressure))/heat_capacity_Cv_B_phase(           T_Tc*temperature_critical_superfluid(pressure),pressure), label='V&W (eq. 3.81) - Winkelmann (eq. 9)')
    plt.xlabel('T/$T_c$')
    plt.ylabel('$Difference$ [J/K/m$^3$]')
    plt.legend()
    plt.show()
