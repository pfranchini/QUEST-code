import numpy as np
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


exec(open("mod_helium3.py").read())

def Width_from_Temperature(Temperature,PressureBar):
    gap = energy_gap_in_low_temp_limit(PressureBar)
    width=np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar)/(2*density*np.pi*d)*np.exp(-gap/(Boltzmann_const*Temperature))
    
    return width


if __name__ == "__main__":

    pressure=0.0
    t_base=100.0e-6 #[K]


    f_base = Width_from_Temperature(1,0)
    noise = np.sqrt(2/(density_of_states(0)*Fermi_velocity(0)*1e-3*1e-6))


    print(noise)
    print(np.sqrt(f_base))
