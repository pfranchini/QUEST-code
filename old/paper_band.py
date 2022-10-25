######################################################
#
# Code to produce plots for the QUEST-DMC paper
#
######################################################

def arguments(argv):

    global pressure
    global d
    arg_help = "{0} -p <pressure [bar]> -d <diameter [m]>".format(argv[0])
    
    try:
        opts, args = getopt.getopt(argv[1:], "hp:d:", ["help", "pressure=", "diameter="])
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
            
def Width_from_Temperature(Temperature,PressureBar,_d):
    """
    Calculates the resonance width [Hz] from the Temperature [K]
    (from Samuli matlab scripts)
    """
    gap = energy_gap_in_low_temp_limit(PressureBar)
    width=np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar)/(2*density*np.pi*_d)*np.exp(-gap/(Boltzmann_const*Temperature))
    
    return width
    
def Temperature_from_Width(Width,PressureBar):
    """
    Calculates the Temperature [K] from the resonance width [Hz]
    (from Samuli matlab scripts)
    """
    
    gap = energy_gap_in_low_temp_limit(PressureBar)
    temperature=-gap/(Boltzmann_const*np.log( Width*2*density*np.pi*d/(np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar))))
    return temperature


############################################


if __name__ == "__main__":
    
    # import Tsepelin code
    execfile("mod_helium3.py")
    
    import sys
    import getopt
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 12})
    
    ####################################################################
    
    pressure = 10  # [bar]
    d = 200e-9;    # [m] vibrating wire
    
    density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)  
    volume = 1e-6      # [m^3] Helium-3 cell
    
    ####################################################################
    
    arguments(sys.argv)
    print('Pressure: ', pressure)
    print('Diameter: ', d)
    
    
    # Dependences on wire diameter, for a certain base temperature and pressure
    t0 = 150e-6   # fix the temperature

    Diameter = np.array([])
    Sensitivity = np.array([])

    for d in np.arange(50e-9, 2000e-9, 100e-9):
        
        W0=Width_from_Temperature(t0,pressure,d) # base width

        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)

        alpha, _ = np.polyfit(DW, DQ, 1)
        Sensitivity = np.append(Sensitivity,1/alpha) # Sensitivity as DeltaDeltaf/DW
        Diameter = np.append(Diameter,d)

    plt.plot(Diameter*1e9,Sensitivity*1e6,label='Sensitivity vs Diameter')
    plt.title('Sensitivity vs Wire diameter ('+str(pressure)+' bar - '+str(t0*1e6)+' $\mu$K)')
    plt.xlabel('diameter [nm]')
    plt.ylabel('$\Delta(\Delta f)$/$\Delta$Q [mHz/keV]')
    #plt.xlim([0, 10e3])   
    plt.savefig('output/Sensitivity_vs_diameter-'+str(int(pressure))+'bar.pdf')
    plt.show()


    # Dependences on T/Tc, for a certain pressure and wire diameter
    d = 100e-9   # fix the wire diameter

    TTc = np.array([])
    Sensitivity = np.array([])

    for ttc in np.arange(0.1, 0.5, 0.01):
        
        W0=Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,d) # base width

        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)

        alpha, _ = np.polyfit(DW, DQ, 1)
        Sensitivity = np.append(Sensitivity,1/alpha) # Sensitivity as DeltaDeltaf/DW
        TTc = np.append(TTc,ttc)


    d = 50e-9   # fix the wire diameter
    TTc = np.array([])
    Sensitivity_d_min = np.array([])
    for ttc in np.arange(0.1, 0.5, 0.01):
        
        W0=Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,d) # base width

        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)

        alpha, _ = np.polyfit(DW, DQ, 1)
        Sensitivity_d_min = np.append(Sensitivity_d_min,1/alpha) # Sensitivity as DeltaDeltaf/DW
        TTc = np.append(TTc,ttc)


    d = 1000e-9   # fix the wire diameter
    TTc = np.array([])
    Sensitivity_d_max = np.array([])
    for ttc in np.arange(0.1, 0.5, 0.01):
        
        W0=Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,d) # base width

        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)

        alpha, _ = np.polyfit(DW, DQ, 1)
        Sensitivity_d_max = np.append(Sensitivity_d_max,1/alpha) # Sensitivity as DeltaDeltaf/DW
        TTc = np.append(TTc,ttc)

    d = 100e-9   # fix the wire diameter
    pressure=0
    TTc = np.array([])
    Sensitivity_p_min = np.array([])
    for ttc in np.arange(0.1, 0.5, 0.01):
        
        W0=Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,d) # base width

        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)

        alpha, _ = np.polyfit(DW, DQ, 1)
        Sensitivity_p_min = np.append(Sensitivity_p_min,1/alpha) # Sensitivity as DeltaDeltaf/DW
        TTc = np.append(TTc,ttc)


    d = 100e-9   # fix the wire diameter
    pressure=30
    TTc = np.array([])
    Sensitivity_p_max = np.array([])
    for ttc in np.arange(0.1, 0.5, 0.01):
        
        W0=Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,d) # base width

        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)

        alpha, _ = np.polyfit(DW, DQ, 1)
        Sensitivity_p_max = np.append(Sensitivity_p_max,1/alpha) # Sensitivity as DeltaDeltaf/DW
        TTc = np.append(TTc,ttc)

    plt.plot(TTc,Sensitivity*1e6,label='Sensitivity vs T/T$_c$')
    plt.fill_between(TTc,Sensitivity_d_min*1e6,Sensitivity_d_max*1e6, hatch="X",facecolor='none', alpha=0.2)
    plt.fill_between(TTc,Sensitivity_p_min*1e6,Sensitivity_p_max*1e6, color='b',alpha=0.2)
    plt.title('Sensitivity vs T/T$_c$ ('+str(pressure)+' bar - '+str(d*1e9)+' nm)')
    plt.xlabel('T/T$_c$')
    plt.ylabel('$\Delta(\Delta f)$/$\Delta$Q [mHz/keV]')
    #ax.xlim([0, 10e3])   
    plt.savefig('output/Sensitivity_vs_TTc-'+str(int(pressure))+'bar.pdf')
    plt.show()

    
    # Dependences on pressure, for a certain temperature and wire diameter
    t0 = 150e-6   # fix the temperature
    d = 100e-9   # fix the wire diameter

    Pressure = np.array([])
    Sensitivity = np.array([])

    for pressure in np.arange(0, 30, 1):
        
        W0=Width_from_Temperature(ttc*temperature_critical_superfluid(pressure),pressure,d) # base width

        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.01):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)

        alpha, _ = np.polyfit(DW, DQ, 1)
        Sensitivity = np.append(Sensitivity,1/alpha) # Sensitivity as DeltaDeltaf/DW
        Pressure = np.append(Pressure,pressure)

    plt.plot(Pressure,Sensitivity*1e6,label='Sensitivity vs Pressure')
    plt.title('Sensitivity vs Pressure ('+str(t0*1e6)+' $\mu$K - '+str(d*1e9)+' nm)')
    plt.xlabel('Pressure [bar]')
    plt.ylabel('$\Delta(\Delta f)$/$\Delta$Q [mHz/keV]')
    #plt.xlim([0, 10e3])   
    plt.savefig('output/Sensitivity_vs_Pressure-.pdf')
    plt.show()




