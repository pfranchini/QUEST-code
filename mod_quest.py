'''
Library to integrate mod_helium3.py for QUEST-DMC calculations
with other derivations reported in the Note https://github.com/pfranchini/QUEST-code/tree/master/note
P.Franchini
'''

def Width_from_Temperature(Temperature,PressureBar,Diameter):
    """
    Calculates the resonance width [Hz] from the Temperature [K]
    (from Samuli matlab scripts)
    eq.10 in the Note
    """
    gap = energy_gap_in_low_temp_limit(PressureBar)
    N_0=density_of_states(PressureBar)/2 # density of quasiparticle states in the normal phase at Fermi energy for one spin component = N_F/2
    width=np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*N_0/(2*density*np.pi*Diameter)*np.exp(-gap/(Boltzmann_const*Temperature))

    return width

def Temperature_from_Width(Width,PressureBar,Diameter):
    """
    Calculates the Temperature [K] from the resonance width [Hz]
    (from Samuli matlab scripts)
    eq.11 in the Note
    """
    gap = energy_gap_in_low_temp_limit(PressureBar)
    N_0=density_of_states(PressureBar)/2 # density of quasiparticle states in the normal phase at Fermi energy for one spin component = N_F/2
    temperature=-gap/(Boltzmann_const*np.log( Width*2*density*np.pi*Diameter/(np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*N_0)
))
    return temperature

def DeltaWidth_from_Energy(E,PressureBar,BaseTemperature,Diameter):
    '''
    Find delta width from the input energy deposition for a certain base temperature
    Find fit line for the Width variation vs Deposited energy for the base temperature
    returns:
     - delta width
     - angular coefficient from the fit
    '''
    
    W0=Width_from_Temperature(BaseTemperature,PressureBar,Diameter)

    DQ = np.array([])  # delta energy [eV]
    DW = np.array([])  # delta width [Hz]

    #for dw in np.arange(0,2.5,0.001):  # Delta(Deltaf)
    for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)  FASTER
        T2= Temperature_from_Width(W0+dw,PressureBar,Diameter)
        T1= Temperature_from_Width(W0,PressureBar,Diameter)
        DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, PressureBar) - heat_capacity_Cv_B_phase_intergral_from_0(T1, PressureBar)) * volume * 6.242e+18) # [eV]
        DW = np.append(DW,dw)

    # Draw the plot
    '''
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
