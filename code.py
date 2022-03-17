
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
            
def Width_from_Temperature(Temperature,PressureBar):
    """
    Calculates the resonance width [Hz] from the Temperature [K]
    (from Samuli matlab scripts)
    """
    gap = energy_gap_in_low_temp_limit(PressureBar)
    width=np.power(Fermi_momentum(PressureBar),2)*Fermi_velocity(PressureBar)*density_of_states(PressureBar)/(2*density*np.pi*d)*np.exp(-gap/(Boltzmann_const*Temperature))
    
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
    
    ####################################################################
    
    pressure = 10  # [bar]
    d = 200e-9;    # [m] vibrating wire
    
    density = 6.05e3;  # [kg/m^3] Niobium-Titanium (NbTi)  
    volume = 1e-6      # [m^3] Helium-3 cell
    
    ####################################################################
    
    arguments(sys.argv)
    print('Pressure: ', pressure)
    print('Diameter: ', d)
    
    
    T = np.array([])   # temperature [K]
    Cv = np.array([])  # heat capacity [J/K/m^3]
    I = np.array([])   # integrated heat capacity (0,T) [J/m^3]
    DQ = np.array([])  # DeltaQ [eV]
    
    for t in np.arange(0.0, 500e-6, 0.000001):
        Cv = np.append(Cv,heat_capacity_Cv_B_phase(t, pressure))
        I  = np.append(I ,heat_capacity_Cv_B_phase_intergral_from_0(t, pressure))
        T  = np.append(T,t)
        
        # Example of deposited energy
        T2 = 150e-6
        deltaT = 0.0005e-3
        DeltaQ = (heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0((T2-deltaT), pressure)) * volume * 6.242e+18 #[eV]
    print("Deposited energy ",DeltaQ," [eV]")
        
        
    #DT = np.array([])  # delta temperature [K]    
    #for dt in np.arange(1e-10, 1e-7, 1e-8):
    #    DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2+dt, pressure) - heat_capacity_Cv_B_phase_intergral_from_0((T2), pressure)) * volume * 6.242e+18) #[eV]#    DT = np.append(DT,dt)
    
    #plt.plot(DQ,DT,label='DQ')
    #plt.title('Temperature variation vs Deposited energy')
    #plt.xlim([0, 100e3])
    #plt.xlabel('$\Delta$Q [eV]')
    #plt.ylabel('$\Delta$T [K]')
    
    
    # Deposited energy for various temperatures at a centain pressure
    T0 = np.array([])  # system temperature [K]
    for t0 in np.arange(40e-6, 210e-6, 20e-6):
        T0 = np.append(T0,t0)
        DQ = np.array([])  # DeltaQ [eV]     
        DT = np.array([])  # delta temperature [K]       
        for dt in np.arange(1e-10, 1e-3, 1e-6):
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(t0+dt, pressure) - heat_capacity_Cv_B_phase_intergral_from_0((t0), pressure)) * volume * 6.242e+18) #[eV] 
            DT = np.append(DT,dt)
        plt.plot(DQ,DT*1e6,label=str(t0*1e6)+' $\mu$K')
                
    plt.title('Temperature variation vs Deposited energy ('+str(pressure)+' bar)')
    plt.xlim([0, 100e3])
    plt.ylim([0, 200])
    plt.xlabel('$\Delta$Q [eV]')
    plt.ylabel('$\Delta$T [$\mu$K]')
    plt.legend()
    plt.savefig('T_vs_DE-'+str(int(pressure))+'bar.pdf')
    plt.show()
    
    # Plots
    #plt.plot(T*1e6,Cv,label='Cv')
    #plt.title('Specific heat capacity vs Temperature')
    #plt.xlim([0, 500])
    #plt.xlabel('T [$\mu$K]')
    #plt.ylabel('$C_v$ [J/K/m$^3$]')
    #plt.savefig('Cv_vs_T-zoom.pdf')
    #plt.show()
    
    plt.plot(T,I,label='I')
    plt.xlim([0, 500e-6])
    plt.title('Integral of heat capacity [0,T] vs Temperature')
    plt.xlabel('T [K]')
    plt.ylabel('I [J/m$^3$]')
    plt.show()
    
    
    T = np.array([])
    W = np.array([])
    
    for w in np.arange(0.0, 200, 0.1):
        
        T = np.append(T,Temperature_from_Width(w,pressure))
        W  = np.append(W,w)
        
    plt.plot(W,T*1e6,label='TvsW')
    plt.title('Temperature vs Resonance width - ('+str(pressure)+' bar - '+str(d*1e9)+' nm)')
    #plt.xlim([0, 500])
    plt.xlabel('$\Delta$f [Hz]')
    plt.ylabel('T [$\mu$K]')
    plt.savefig('T_vs_W.pdf')
    plt.show()

    
    T = np.array([])
    W = np.array([])
    
    for t in np.arange(50e-6, 300e-6, 1e-6):
        
        W = np.append(W,Width_from_Temperature(t,pressure))
        T  = np.append(T,t)
        
    plt.plot(T*1e6,W,label='WvsT')
    plt.title('Resonance width vs Temperature - ('+str(pressure)+' bar - '+str(d*1e9)+' nm)')
    #plt.xlim([0, 500])
    plt.xlabel('T [$\mu$K]')
    plt.ylabel('$\Delta$f [Hz]')
    plt.savefig('W_vs_T.pdf')
    plt.show()
    
    
    # Width variation vs Deposited energy for multiple starting temperatures
    for t0 in np.arange(50e-6, 310e-6, 50e-6):
        
        W0=Width_from_Temperature(t0,pressure)
        print(t0,W0)
        
        DQ = np.array([])  # delta energy [eV]
        DW = np.array([])  # delta width [Hz]    
    
        for dw in np.arange(0,2.5,0.001):  # Delta(Deltaf)
            T2= Temperature_from_Width(W0+dw,pressure)
            T1= Temperature_from_Width(W0,pressure)
            #print(T1,T2,T2-T1)
            DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
            DW = np.append(DW,dw)
            #print(dw,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18)

        # Draw a line for each initial temperature
        plt.plot(DQ,DW,label=str(t0*1e6)+' $\mu$K')

        # Fit each line
        m, q = np.polyfit(DW, DQ, 1)
        print("m: ",m)
        
    plt.title('Width variation vs Deposited energy ('+str(pressure)+' bar - '+str(d*1e9)+' nm)')
    plt.xlim([0, 100e3])
    plt.ylim([0, 0.3])
    plt.xlabel('$\Delta$Q [eV]')
    plt.ylabel('$\Delta(\Delta f)$ [Hz]')
    plt.legend()
    plt.savefig('DeltaDeltaW_vs_DE-'+str(int(pressure))+'bar.pdf')
    plt.show()



    # Bolometric calibration (Winkelmann)
    pressure=0
    d=4.5e-06
    volume=3*0.1413*1e-6
    # Width variation vs Deposited energy 
    W0=1.7 # Base width as in the paper
    print("Base temperature: ",str(Temperature_from_Width(W0, pressure)))
    DQ = np.array([])  # delta energy [eV]
    DW = np.array([])  # delta width [Hz]    
    
    for dw in np.arange(0,0.3,0.01):
        T2= Temperature_from_Width(W0+dw,pressure)
        T1= Temperature_from_Width(W0,pressure)
        #print(T1,T2,T2-T1)
        DQ = np.append(DQ,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18) #[eV]#
        DW = np.append(DW,dw)
        #print(dw,(heat_capacity_Cv_B_phase_intergral_from_0(T2, pressure) - heat_capacity_Cv_B_phase_intergral_from_0(T1, pressure)) * volume * 6.242e+18)
                
    plt.plot(DQ/1e3,DW*1e3,label='DQvsDW')
    plt.title('Width variation vs Deposited energy (Winkelmann)')
    plt.xlim([0, 900])
    plt.xlabel('$\Delta$Q [KeV]')
    plt.ylabel('$\Delta(\Delta f)$ [mHz]')
    plt.savefig('Winkelmann.pdf')
    plt.show()
