# -*- coding: utf-8 -*-
"""
Created on Sat June 03 2017
library should take in arguments as float, list, numpy arrays or Pandas
@author: Tsepelin V
"""
import numpy as np
#import mod_oscillator as osc
#import mod_stokes as stokes
#import mod_fileio as fio
from scipy.special import erfc

##############################################################################
###     Constants
##############################################################################


Boltzmann_const = 1.38064852e-23 #[J/K]
Plank_const = 6.626176e-34
Plankbar_const = Plank_const / (2.0 * np.pi)
Avogadro_const = 6.022140857e+23
ElectronV = 1.6021766208e-19 #[C]
atomic_mass_unit = 1.660539040e-27 #[kg]

molar_mass  = 3.0160293e-3 #[kg]
atomic_mass = 3.0160293 * atomic_mass_unit
kappa = Plank_const /(2 * molar_mass / Avogadro_const)
gap_coeff = 1.76 # BCS theory energy gap value

#coefficients used in "mod_stokes" for vibrating wire thermometry
meanfreepath_adjustment = 1.9
ballistic_switch = 1.0

#superfluid coefficents for thermometry abd BBR calibration
gamma_coefficient = 0.28 # PRB 57 (1998) 14381, Bauerle et.al. Grenoble + Lancaster

yoshida_scaledT_poly = [-0.0611, 0.1396, - 0.0216, -0.0075, 0.9074]


##############################################################################
###     Normal liquid parameters
##############################################################################

def molarvolume(PressureBar):
    """For a given pressure [bar] function returns the helium-3 molar volume [m^3 per mol].

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    molarvolume_poly = [-0.91253577e-6, 0.94759780e-4,  -0.38859562e-2, 0.83421417e-1, -0.11803474e1, 36.837231]

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return 1.0e-6 * np.polyval(molarvolume_poly, PressureBar)

def he3_density(PressureBar):
    """For a given pressure[bar] function returns the helium-3 density [kg m^3].

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return molar_mass / molarvolume(PressureBar)


def density_SVP(TemperatureK):
    """For a given temperature[K] function returns the helium-3 density at SVP [kg m^3].

    based on Willks?
    ; http://dx.doi.org/
   ''
   """

    density_poly = [0.0604, - 0.1175,  0.8375, - 0.7880, 36.9087]

    # change float into numpy array for function len() to work
    try:
        len(TemperatureK)
    except:
        TemperatureK = np.array([TemperatureK])

    rho = 1.0e6* molar_mass / np.polyval(density_poly, TemperatureK)

    return rho

def sound_velocity_first(TemperatureK):
    """
    Calculates the first sound velocity of helium-3 for a given temperature

    :param TemperatureK: The temperature in Kelvin
    :return: The first sound velocity in m/s

	based on Laquer H.L. et al., Phys. Rev. 113, 417 (1959)
    """
    first_sound_poly = [-0.00176, 0.0, 0.0, 0.0, 0.0, -0.130, -5.98, 0.0, 183.9]

    try:
        len(TemperatureK)
    except TypeError as e:
        TemperatureK = np.array([TemperatureK])

    first_sound = np.zeros(len(TemperatureK))

    for temp_index in range(len(TemperatureK)):
        first_sound[temp_index] = np.polyval(first_sound_poly, TemperatureK)

    return first_sound


def heatcapacity_gamma(PressureBar):
    """For a given pressure[bar] function returns the value of heat capacity gamma .

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    gamma_poly = [-0.53785385e-6,0.46153498e-4,-0.14738303e-2,0.69575243e-1,2.7840464]

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return np.polyval(gamma_poly, PressureBar)

def mass_effective(PressureBar):
    """For a given pressure[bar] function returns the value of effective mass.

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return 3.06413e-27 * 1.0e-2 * heatcapacity_gamma(PressureBar) / (molarvolume(PressureBar) * Fermi_momentum(PressureBar))


def mass_effective_polynomial(PressureBar):
    """For a given pressure[bar] function returns the value of effective mass.

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    masseff_poly = [ 1.75334e-8, -2.4401e-6, 1.33671e-4, -0.00366, 0.13133, 2.79972]

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return np.polyval(masseff_poly, PressureBar)


def viscosity(TemperatureK, PressureBar = 0.0):
    """For a given pressure[Bar] and temperature[K] array returns the helium-3 viscosity [Pa s].

    Function is valid between 5?mK and 17 mK
    based on D. C. Carless,t H. E. Hall, and J. R. Hook
    Journal of Low Temperature Physics 50, 584 (1983); http://dx.doi.org/
    'Vibrating Wire Measurements in Liquid 3He. I. The Normal State'
    1.0e-1 / (coeff_A*(1.0e3*TemperatureK)**2 + coeff_B)
    """

    if PressureBar < 1.28:
        tempA = 0.38 - 0.007 * (1.28 - PressureBar) / 1.18
        tempB = 0.06 - 0.02 * (1.28 - PressureBar) / 1.18
    elif PressureBar < 4.65:
        tempA = 0.424 - 0.044 * (4.65 - PressureBar) / 3.37
        tempB = 0.19 - 0.13 * (4.65 - PressureBar) / 3.37
    elif PressureBar < 9.89:
        tempA = 0.495 - 0.071 * (9.89 - PressureBar) / 5.24
        tempB = 0.43 - 0.24 * (9.89 - PressureBar) / 5.24
    elif PressureBar < 19.89:
        tempA = 0.603 - 0.108 * (19.89 - PressureBar) / 10
        tempB = 0.94 - 0.56 * (19.89 - PressureBar) / 10
    else:
        tempA = 0.603 + 0.107 * (PressureBar - 19.89) / 9.45
        tempB = 0.94 + 0.56 * (PressureBar - 19.89) / 9.45

    coeffA = tempA * 10 * 1.12e3 * 1.12e3
    coeffB = tempB * 10

    # change float into numpy array for function len() to work
    try:
        len(TemperatureK)
    except:
        TemperatureK = np.array([TemperatureK])

    return 1.0 / (coeffA * TemperatureK**2 + coeffB)  # Pa s


def viscosity_SVP(TemperatureK):
    """For a given temperature[K] array returns the helium-3 viscosity at SVP[Pa s].

    Function is valid between ~6mK and 1.8K
    based on D.I. Bradley et al.
    J Low Temp Phys 171, 750 (2013) ; http://dx.doi.org/10.1007/s10909-012-0804-3
    'Thermometry in Normal Liquid 3He Using a Quartz Tuning Fork Viscometer'
    """
    viscosity_poly = [-6.92457,22.39929,-24.7364,28.22503,4.73534,2.02717]

    # change float into numpy array for function len() to work
    try:
        len(TemperatureK)
    except:
        TemperatureK = np.array([TemperatureK])

    return 1.0e-7 * np.polyval(viscosity_poly,TemperatureK) / np.power(TemperatureK, 2.0) # Pa s

def viscosity_SVP_BHT(TemperatureK):
    """For a given temperature[K] array returns the helium-3 viscosity [Pa s] at SVP.

    Function is valid between 50mK and 3.0K
    based on M. A. BLACK, H. E. HALL and K. THOMPSON
    J. Phys. C: Solid St. Phys., 4, 129 (1971); http://dx.doi.org/
    'The viscosity of liquid helium 3'
    """

    # change float into numpy array for function len() to work
    try:
        len(TemperatureK)
    except:
        TemperatureK = np.array([TemperatureK])

    return 1.0e-7 * (2.21 / TemperatureK**2 + 26.3 / TemperatureK**(1/3) ) # Pa s

def viscosity_SVP_CHH(TemperatureK):
    """For a given temperature[K] array returns the helium-3 viscosity [Pa s] at SVP.

    Function is valid between 5?mK and 17 mK
    based on D. C. Carless,t H. E. Hall, and J. R. Hook
    Journal of Low Temperature Physics 50, 584 (1983); http://dx.doi.org/
    'Vibrating Wire Measurements in Liquid 3He. I. The Normal State'
    """

    # change float into numpy array for function len() to work
    try:
        len(TemperatureK)
    except:
        TemperatureK = np.array([TemperatureK])

    return 1.0e-1 / (0.468*(1.0e3*TemperatureK)**2 + 0.0383)  # Pa s

def viscosity_Tony(TemperatureK):
    if type(TemperatureK) == float:
        TemperatureK = np.array([TemperatureK])
    eta = []
    for temp in TemperatureK.tolist():
        e = 0.277e-7 / temp**2 + 3.4e-7 / temp
        if temp < 0.0165:
            e = 0.305e-7 / temp**2 + 1.35e-7 / temp + 2.2e-6

        if temp > 0.068:
            e = 0.29e-7 / temp**2 + 1.65e-7 / temp + 2.3e-6
        eta.append(e)

    return np.array(eta)


def viscosity_over_meanfreepath(PressureBar):
    """For a given temperature[K] and Concentration returns the He3-He4 mixture mean free path [???].

    Function is valid between ~?mK and ?K
    based on D.I. Bradley et al.
    """

    return 0.2 * (Avogadro_const/molarvolume(PressureBar)) ** (4.0/3.0) * (3.0*9.8696) ** (1.0/3.0) * Plankbar_const

def normalfluid_sliplength(TemperatureK, PressureBar):
    """
    For a given temperatre [k] and pressure [bar] returns the slip length of He-3.
    Equation from Dieter Vollhardt, Peter Wolfle "The Superfluid Phases of Helium-3" pp 483
    slip length = 0.579 * mean free path

    :param TemperatureK: Temperature in Kelvin
    :param PressureBar: Pressure in Bar
    :return: The slip length in meters
    """
    return 0.579 * meanfreepath(TemperatureK, PressureBar)

def meanfreepath(TemperatureK, PressureBar):
    """For a given temperature[K] and Concentration returns the He3-He4 mixture mean free path [???].

    Function is valid between ~?mK and ?K
    based on D.I. Bradley et al.
    """

    return viscosity(TemperatureK, PressureBar) / viscosity_over_meanfreepath(PressureBar)

def meanfreepath_SVP(TemperatureK, PressureBar = 0.0):
    """For a given temperature[K] and Concentration returns the He3-He4 mixture mean free path [???].

    Function is valid between ~?mK and ?K
    based on D.I. Bradley et al.
    """

    return viscosity_SVP(TemperatureK) / viscosity_over_meanfreepath(PressureBar)

def Fermi_momentum(PressureBar):
    """For a given pressure[bar] function returns the Fermi momentum of helium-3 [kg m s-1].

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return 2.7551e-26 * molarvolume(PressureBar)**(-1.0/3.0)

def Fermi_velocity(PressureBar):
    """For a given pressure[bar] function returns the Fermi velocity of helium-3 [m s-1].

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return Fermi_momentum(PressureBar) / (mass_effective(PressureBar) * atomic_mass)

def Fermi_parameterF1(PressureBar):
    """For a given pressure[bar] function returns the value of F1 Fermi parameter.

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return 3.0 * (mass_effective(PressureBar) - 1)

def density_of_states(PressureBar):
    """For a given pressure[Bar] returns the density of states at Fermi energy

    Function is valid between 0 and 35bars
    based on Enns Eq.(3.28) chapter on Landau Fermi Theory
    """

    return 2.0**3 * np.pi * atomic_mass * mass_effective(PressureBar) * Fermi_momentum(PressureBar) / Plank_const ** 3.0

def K0(PressureBar):
    """For a given pressure[Bar] returns the K0[???].

    Function is valid between 0 and 35bars
    based on D.I. Bradley et al.
    """

    return 8.32 / 3.0 * heatcapacity_gamma(PressureBar) * Fermi_velocity(PressureBar) / molarvolume(PressureBar)

def Kcorr(PressureBar):
    """For a given pressure[Bar] returns the Kcorr[???].

    Function is valid between 0 and 35bars
    based on D.I. Bradley et al.
    """

    return np.pi * 0.95e-3**2.0 * Boltzmann_const**2.0 / 2.0 * density_of_states(PressureBar)* Fermi_velocity(PressureBar)


def viscous_penetration_depth(FrequencyHz, TemperatureK, PressureBar):
    """
-    Calculates the viscous penetration depth for a given oscillator with frequency(Hz) in He3 fluid.
-
-    Uses delta = sqrt(eta/(density_nf * pi * omega)
-
-    FrequencyHz: frequency of oscillator
-    TemperatureK: temperature of helium mixture at SVP
-    """

    deltaM = osc.viscous_penetration_depth(2.0 * FrequencyHz, viscosity(TemperatureK, PressureBar), density(PressureBar))
    #double frequency is used to be consistent with an old library

    return deltaM


def viscous_penetration_depth_SVP(FrequencyHz, TemperatureK):
    """
-    Calculates the viscous penetration depth for a given oscillator with frequency(Hz) in He3He4 mixture.
-
-    Uses delta = sqrt(eta/(density_nf * pi * omega)
-
-    FrequencyHz: frequency of oscillator
-    TemperatureK: temperature of helium mixture at SVP
-    """

    deltaM = osc.viscous_penetration_depth(2.0 * FrequencyHz, viscosity_SVP(TemperatureK), density_SVP(TemperatureK))
    #double frequency is used to be consistent with an old library

    return deltaM

def __get_charatersictic_size(ObjOscillator):
    if type(ObjOscillator) == osc.vibrating_loop:
        length = ObjOscillator._diameter
    elif type(ObjOscillator) == osc.tuning_fork:
        length = ObjOscillator.width()
    else:
        raise NotImplementedError('Object oscillator "%s" has no implemented characteristic size' %
                                  ObjOscillator.name())
    return length / ObjOscillator.get_stokes_parameter()


def resonance_frequency_shift_obj(ObjOscillator, TemperatureK, PressureBar):
    """
-    Calculates frequency shift of an oscillator in normal He3 for a given oscillator.
-
-    Made from the old fortran code by Bradley, parameters are tuned for the Vibrating wires
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP, the module will not work correctly above ~100mK
-    """

    osc_radius = __get_charatersictic_size(ObjOscillator)
    osc_frequency = ObjOscillator.frequency_vacuum
    osc_density = ObjOscillator.density()

    res_shift = stokes.resonance_frequency_shift(osc_frequency, osc_radius, osc_density,
                                viscous_penetration_depth(osc_frequency, TemperatureK, PressureBar),
                                meanfreepath(TemperatureK, PressureBar),
                                density(PressureBar), meanfreepath_adjustment, ballistic_switch)

    return res_shift

def resonance_frequency_shift_obj_SVP(ObjOscillator, TemperatureK, PressureBar = 0.0):
    """
-    Calculates frequency shift of an oscillator in normal He3 for a given oscillator.
-
-    Made from the old fortran code by Bradley, parameters are tuned for the Vibrating wires
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP, the module will not work correctly above ~100mK
-    """

    osc_radius = __get_charatersictic_size(ObjOscillator)
    osc_frequency = ObjOscillator.frequency_vacuum
    osc_density = ObjOscillator.density()

    res_shift = stokes.resonance_frequency_shift(osc_frequency, osc_radius, osc_density,
                                viscous_penetration_depth_SVP(osc_frequency, TemperatureK),
                                meanfreepath_SVP(TemperatureK),
                                density(PressureBar), meanfreepath_adjustment, ballistic_switch)

    return res_shift


def resonance_width_obj(ObjOscillator, TemperatureK, PressureBar):
    """
-    Calculates frequency width (damping) of an oscillator in dilute mixture for a given oscillator.
-
-    Made from the old fortran code by Bradley, parameters are tuned for the Vibrating wires
    Much more sophistificated than a standard Stokes model
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP, the module will not work correctly above ~100mK
-    """

    osc_radius = __get_charatersictic_size(ObjOscillator)
    osc_frequency = ObjOscillator.frequency_vacuum
    osc_density = ObjOscillator.density()

    res_width = stokes.resonance_width(osc_frequency, osc_radius, osc_density,
                                viscous_penetration_depth(osc_frequency, TemperatureK, PressureBar),
                                meanfreepath(TemperatureK, PressureBar),
                                density(PressureBar), meanfreepath_adjustment, ballistic_switch)

    return res_width

def resonance_width_obj_SVP(ObjOscillator, TemperatureK, PressureBar = 0.0):
    """
-    Calculates frequency width (damping) of an oscillator in dilute mixture for a given oscillator.
-
-    Made from the old fortran code by Bradley, parameters are tuned for the Vibrating wires
    Much more sophistificated than a standard Stokes model
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP, the module will not work correctly above ~100mK
-    """

    osc_radius = __get_charatersictic_size(ObjOscillator)
    osc_frequency = ObjOscillator.frequency_vacuum
    osc_density = ObjOscillator.density()

    res_width = stokes.resonance_width(osc_frequency, osc_radius, osc_density,
                                viscous_penetration_depth_SVP(osc_frequency, TemperatureK),
                                meanfreepath_SVP(TemperatureK),
                                density(PressureBar), meanfreepath_adjustment, ballistic_switch)

    return res_width



def temperature_from_width_obj(ObjOscillator, currentWidthHz, PressureBar):
    """
-    Calculates the temperature of 3He normal fluid from the current width for a given oscillator.
-
-    Made from old fortran code by Bradley, parameters are tuned for the Vibrating wires
    Much more sophisticated than a standard Stokes model
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP
-    """

    # define required temperature range, the module will not work correctly above ~100mK
    TemperatureK = np.linspace(temperature_critical_superfluid(PressureBar), 0.2, num=10000) # generate temperature array
    # Calculate resonance width for all temperatures
    res_width = resonance_width_obj(ObjOscillator, TemperatureK, PressureBar)

    # change float into numpy array for function len() to work
    try:
        len(currentWidthHz)
    except:
        currentWidthHz = np.array([currentWidthHz])
    currentWidthHz = np.array(currentWidthHz)

    currentTemperatureK = np.tile(np.NaN,len(currentWidthHz)) # array of NaN

    for tempindex in range(len(currentTemperatureK)):
        # Find the closest width an hence the temperature
        minindex = np.argmin(np.abs(res_width - np.tile(currentWidthHz[tempindex], len(TemperatureK))))
        currentTemperatureK[tempindex] = TemperatureK[minindex]

    return currentTemperatureK

def temperature_from_width_obj_SVP(ObjOscillator, currentWidthHz, PressureBar = 0.0):
    """
-    Calculates the temperature of 3He normal fluid from the current width for a given oscillator.
-
-    Made from old fortran code by Bradley, parameters are tuned for the Vibrating wires
    Much more sophisticated than a standard Stokes model
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP
-    """

    Max_valid_temperature = 1.6 #[K]

    # define required temperature range, the module will not work correctly above ~100mK
    TemperatureK = np.linspace(temperature_critical_superfluid(PressureBar),Max_valid_temperature,num=60000) # generate temperature array
    # Calculate resonance width for all temperatures
    res_width = resonance_width_obj(ObjOscillator, TemperatureK, PressureBar)

    # change float into numpy array for function len() to work
    try:
        len(currentWidthHz)
    except:
        currentWidthHz = np.array([currentWidthHz])
    currentWidthHz = np.array(currentWidthHz)

    currentTemperatureK = np.tile(np.NaN,len(currentWidthHz)) # array of NaN

    for tempindex in range(len(currentTemperatureK)):
        # Find the closest width an hence the temperature
        minindex = np.argmin(np.abs(res_width - np.tile(currentWidthHz[tempindex], len(TemperatureK))))
        currentTemperatureK[tempindex] = TemperatureK[minindex]

    return currentTemperatureK



##############################################################################
###     Superfluid liquid parameters
##############################################################################

def energy_gap_in_low_temp_limit(PressureBar):
    return gap_coeff * temperature_critical_superfluid(PressureBar) * Boltzmann_const

def real_squashing_energy_in_low_temp_limit(PressureBar):
    return np.sqrt(8./5.) * energy_gap_in_low_temp_limit(PressureBar)

def imaginary_squashing_energy_in_low_temp_limit(PressureBar):
    return np.sqrt(12./5.) * energy_gap_in_low_temp_limit(PressureBar)

def pair_breaking_energy_in_low_temp_limit(PressureBar):
    return 2 * energy_gap_in_low_temp_limit(PressureBar)

def pair_breaking_superfluid_velocity_in_low_temp_limit(PressureBar, energy_gap_suppression=1.0):
    return energy_gap_suppression * energy_gap_in_low_temp_limit(PressureBar) / (Fermi_momentum(PressureBar))

def real_squashing_superfluid_velocity_in_low_temp_limit(PressureBar, energy_gap_suppression=1.0):
    return (energy_gap_suppression * np.sqrt(5. / 8.) * energy_gap_in_low_temp_limit(PressureBar) /
            Fermi_momentum(PressureBar))

def imaginary_squashing_superfluid_velocity_in_low_temp_limit(PressureBar, energy_gap_suppression=1.0):
    return (energy_gap_suppression * np.sqrt(5. / 12.) * energy_gap_in_low_temp_limit(PressureBar) /
            Fermi_momentum(PressureBar))

def Landau_velocity(PressureBar):
    """For a given pressure[bar] function returns the Landau velocity of helium-3 [m s-1].

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return energy_gap_in_low_temp_limit(PressureBar) / Fermi_momentum(PressureBar)


def density_superfluid(TemperatureK, PressureBar):
    """For a given pressure[bar] function returns the effective helium-3 density in superfluid [kg m^3].

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

#    print(mass_effective(PressureBar) * yoshidafunction_0(TemperatureK / temperature_critical_superfluid(PressureBar)) / \
#        (1.0 + Fermi_parameterF1(PressureBar) * yoshidafunction_0(TemperatureK / temperature_critical_superfluid(PressureBar)) / 3.0))
    density_effective = density(PressureBar) * mass_effective(PressureBar) * yoshidafunction_0(TemperatureK / temperature_critical_superfluid(PressureBar)) / \
        (1.0 + Fermi_parameterF1(PressureBar) * yoshidafunction_0(TemperatureK / temperature_critical_superfluid(PressureBar)) / 3.0)

    return density_effective


def temperature_critical_superfluid(PressureBar):
    """For a given pressure[bar] function returns the transition temperature for superfluidity [K].

    based on Greywall
    ; http://dx.doi.org/
   ''
   """

    CritTempPoly = [0.53010918e-7,-0.57248644e-5,0.25685169e-3,-0.69302185e-2,0.13867188,0.92938375]

    # change float into numpy array for function len() to work
    try:
        len(PressureBar)
    except:
        PressureBar = np.array([PressureBar])

    return 1.0e-3 * np.polyval(CritTempPoly, PressureBar)

def coherence_length(PressureBar):
    """For a given pressure[Bar] returns the coherence length [m].

    Function is valid between 0 and 35bars
    based on D.I. Bradley et al.
    """

    return Plankbar_const * Fermi_velocity(PressureBar) / (2.0 * np.pi* Boltzmann_const * temperature_critical_superfluid(PressureBar))

def yoshidafunction_0(reduced_T):
    """
    % yosida functions Y0
    % CHH expression: t is reduced temperature (T/Tc)
    % ts is their scaled temperature
    """
    # change float into numpy array for function len() to work
    try:
        len(reduced_T)
    except:
        reduced_T = np.array([reduced_T])

    Tbelow0p94 = [-1.367, 4.625, -0.88, 3.454]
    Tabove0p94 = [1.985, -0.985]

    scaledT = np.polyval(yoshida_scaledT_poly, reduced_T) * reduced_T

    yoshidabelow0p94 = np.exp(-1.76388 / scaledT) * np.polyval(Tbelow0p94, scaledT) / np.sqrt(scaledT)
    yoshidaabove0p94 = np.polyval(Tabove0p94, scaledT)

    return (reduced_T >= np.tile(0.94,len(reduced_T))) * yoshidaabove0p94 + (reduced_T < np.tile(0.94,len(reduced_T))) * yoshidabelow0p94

def yoshidafunction_5(reduced_T):
    """
    % yosida functions Y5
    % CHH expression: t is reduced temperature (T/Tc)
    % ts is their scaled temperature
    """

    # change float into numpy array for function len() to work
    try:
        len(reduced_T)
    except:
        reduced_T = np.array([reduced_T])

    Tbelow0p80 = [0.392, -1.425, 1.1958, 0.10177]

    scaledT = np.polyval(yoshida_scaledT_poly, reduced_T) * reduced_T

    yoshidabelow0p80 = np.exp(1.76388 / scaledT) * np.polyval(Tbelow0p80, scaledT) / np.sqrt(scaledT)
    yoshidaabove0p80 = np.exp(1.76388 / scaledT) * (0.19847 + 0.335 * np.sqrt(1-scaledT)) / np.sqrt(scaledT)

    return (reduced_T >= np.tile(0.80,len(reduced_T))) * yoshidaabove0p80 + (reduced_T < np.tile(0.80,len(reduced_T))) * yoshidabelow0p80

def yoshidafunction_6(reduced_T):
    """
    % yosida functions Y5
    % CHH expression: t is reduced temperature (T/Tc)
    % ts is their scaled temperature
    """

    # change float into numpy array for function len() to work
    try:
        len(reduced_T)
    except:
        reduced_T = np.array([reduced_T])

    Tbelow0p90 = [4.1, -2.117, 0.4467, 2.402]
    Tabove0p90 = [-7.5, 13.275, -4.517]

    scaledT = np.polyval(yoshida_scaledT_poly, reduced_T) * reduced_T

    yoshidabelow0p90 = np.exp(-1.76388 / scaledT) * np.polyval(Tbelow0p90, scaledT)
    yoshidaabove0p90 = 1 - np.sqrt(1 - scaledT) * np.polyval(Tabove0p90, scaledT)

    return (reduced_T >= np.tile(0.90,len(reduced_T))) * yoshidaabove0p90 + (reduced_T < np.tile(0.90,len(reduced_T))) * yoshidabelow0p90

def heat_capacity_Cv_B_phase(TemperatureK, PressureBar):
    """
    Calculates the heat capacity under constant volume. Uses equation from Vollhardt and Wolfle: The Superfluid Phases
    of Helium 3 p84.
    :param TemperatureK: The temperature in Kelvin
    :param PressureBar: The pressure in bar
    :return: The heat capacity in J / m ^3 / K
    """
    gap = energy_gap_in_low_temp_limit(PressureBar)
    cv = np.sqrt(2 * np.pi) * Boltzmann_const * density_of_states(PressureBar) * gap
    cv = cv * np.power(gap / (Boltzmann_const * TemperatureK), 3./2.) * np.exp(-1 * gap / (Boltzmann_const * TemperatureK))
    return cv # J / m3 / K

def heat_capacity_Cv_B_phase_intergral_from_0(TemperatureK, PressureBar):
    """
    Returns the definite integral of the heat capacity from absolute zero with respect to temperature. Uses equation
    from Vollhardt and Wolfle: The Superfluid Phases of Helium 3 p84.
    :param TemperatureK: The temperature in Kelvin
    :param PressureBar: The pressure in bar
    :return: The energy per unit volume in J / m^3
    """
    gap = energy_gap_in_low_temp_limit(PressureBar)
    cv = np.sqrt(2 * np.pi) * Boltzmann_const * density_of_states(PressureBar) * gap
    intergral = cv *  np.sqrt(np.pi) * gap / Boltzmann_const * erfc(np.sqrt(gap / (Boltzmann_const * TemperatureK)))
    return intergral #  J / m3

def viscosity_superfluid(TemperatureK, PressureBar = 0.0):
    """For a given pressure[Bar] and temperature[K] array returns the effective viscosity in superfluid helium-3 [Pa s].

    based on tony's Shaun's calculations

    only works for T < superfluidTc

    a plot of CHH  reduced viscosity gives the following reasonable
    form. Basically, LOG( redvis) vz sqrt(1-T)^.5 reasonably slow

    """
    # change float into numpy array for function len() to work
    try:
        len(TemperatureK)
    except:
        TemperatureK = np.array([TemperatureK])

    redT = TemperatureK / temperature_critical_superfluid(PressureBar)

    redviscoeff_bl0p6 = 0.11
    redviscoeff_bl0p7 = 10.0**(-0.9586-0.8562*(np.sqrt(1-redT)-np.sqrt(0.4)))
    redviscoeff_bl0p8 = 10.0**(-0.8861-0.6183*(np.sqrt(1-redT)-np.sqrt(0.3)))
    redviscoeff_bl0p9 = 10.0**(-0.8239-1.4172*(np.sqrt(1-redT)-np.sqrt(0.2)))
    redviscoeff_bl0p95 = 10.0**(-0.6383-1.7352*(np.sqrt(1-redT)-np.sqrt(0.1)))
    redviscoeff_bl0p975 = 10.0**(-0.4776-1.6177*(np.sqrt(1-redT)-np.sqrt(0.05)))
    redviscoeff_bl1p0 = 10.0**(-0.3716-2.3503*(np.sqrt(1-redT)-np.sqrt(0.025)))

    redvis = (redT <= np.tile(0.6,len(redT))) * redviscoeff_bl0p6 + \
            (redT >= np.tile(0.6,len(redT))) * (redT < np.tile(0.7,len(redT))) * redviscoeff_bl0p7 + \
            (redT >= np.tile(0.7,len(redT))) * (redT < np.tile(0.8,len(redT))) * redviscoeff_bl0p8 + \
            (redT >= np.tile(0.8,len(redT))) * (redT < np.tile(0.9,len(redT))) * redviscoeff_bl0p9 + \
            (redT >= np.tile(0.9,len(redT))) * (redT < np.tile(0.95,len(redT))) * redviscoeff_bl0p95 + \
            (redT >= np.tile(0.95,len(redT))) * (redT < np.tile(0.975,len(redT))) * redviscoeff_bl0p975 + \
            (redT >= np.tile(0.975,len(redT))) * (redT < np.tile(1.0,len(redT))) * redviscoeff_bl1p0 + \
            (redT >= np.tile(1.0,len(redT))) * 1.0

    return redvis * viscosity(temperature_critical_superfluid(PressureBar), PressureBar)


def slip_length_superfluid(TemperatureK, PressureBar):
    """
    % zeta functions: effective slip length
    % CHH expression: t is reduced temperature (T/Tc)
    % ts is their scaled temperature
    """

    return 0.5*yoshidafunction_5(TemperatureK/temperature_critical_superfluid(PressureBar)) * viscosity_superfluid(TemperatureK) /\
           viscosity_over_meanfreepath(PressureBar)

def group_velocity_superfluid(TemperatureK, PressureBar):
    """
    Returns the approximate group velocity for helium-3 superfluid B
    :param TemperatureK: Temperature in Kelvin
    :param PressureBar: Pressure in bar
    :return: Approximate group velocity in meters per second
    """
    return Fermi_velocity(PressureBar) * np.sqrt(TemperatureK /
                                                 (gap_coeff * temperature_critical_superfluid(PressureBar)))

def meanfreepath_superfluid(TemperatureK, PressureBar):
    """For a given temperature[K] and pressure returns mean free path in superfluid He3 [???].

    Function is valid between ~?mK and ?K
    based on D.I. Bradley et al.

    """

    return viscosity_superfluid(TemperatureK, PressureBar) / \
           (viscosity_over_meanfreepath(PressureBar) * yoshidafunction_6(TemperatureK/temperature_critical_superfluid(PressureBar)))

def meanfreepath_adjustment_superfluid(TemperatureK, PressureBar):
    """For a given temperature[K] and pressure returns mean free path adjustment (fudge factor alpha) in superfluid He3 [???].

    Function is valid between ~?mK and ?K
    based on D.I. Bradley et al.

    """

    return 1.156 * meanfreepath_adjustment / (yoshidafunction_5(TemperatureK/temperature_critical_superfluid(PressureBar))*yoshidafunction_6(TemperatureK/temperature_critical_superfluid(PressureBar)))


def viscous_penetration_depth_superfluid(FrequencyHz, TemperatureK, PressureBar):
    """
-    Calculates the viscous penetration depth for a given oscillator with frequency(Hz) in He3 fluid.
-
-    Uses delta = sqrt(eta/(density_nf * pi * omega)
-
-    FrequencyHz: frequency of oscillator
-    TemperatureK: temperature of helium mixture at SVP
-    """

    deltaM = osc.viscous_penetration_depth(2.0 * FrequencyHz, viscosity_superfluid(TemperatureK, PressureBar),
                                           density_superfluid(TemperatureK, PressureBar))
    #double frequency is used to be consistent with an old library

    return deltaM


def resonance_frequency_shift_superfluid_obj(ObjOscillator, TemperatureK, PressureBar):
    """
-    Calculates frequency shift of an oscillator in normal He3 for a given oscillator.
-
-    Made from the old fortran code by Bradley, parameters are tuned for the Vibrating wires
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP, the module will not work correctly above ~100mK
-    """

    osc_radius = ObjOscillator._diameter / 2.0
    osc_frequency = ObjOscillator.frequency_vacuum
    osc_density = ObjOscillator.density()

    res_shift = stokes.resonance_frequency_shift_superfluidHe3B(osc_frequency, osc_radius, osc_density,
                                viscous_penetration_depth_superfluid(osc_frequency, TemperatureK, PressureBar),
                                meanfreepath_superfluid(TemperatureK, PressureBar),
                                slip_length_superfluid(TemperatureK, PressureBar),
                                density_superfluid(TemperatureK, PressureBar),
                                meanfreepath_adjustment_superfluid(TemperatureK, PressureBar), ballistic_switch)

    return res_shift


def resonance_width_superfluid_obj(ObjOscillator, TemperatureK, PressureBar):
    """
-    Calculates frequency width (damping) of an oscillator in dilute mixture for a given oscillator.
-
-    Made from the old fortran code by Bradley, parameters are tuned for the Vibrating wires
    Much more sophistificated than a standard Stokes model
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP, the module will not work correctly above ~100mK
-    """

    osc_radius = ObjOscillator._diameter / 2.0
    osc_frequency = ObjOscillator.frequency_vacuum
    osc_density = ObjOscillator.density()

    res_width = stokes.resonance_width_superfluidHe3B(osc_frequency, osc_radius, osc_density,
                                viscous_penetration_depth_superfluid(osc_frequency, TemperatureK, PressureBar),
                                meanfreepath_superfluid(TemperatureK, PressureBar),
                                slip_length_superfluid(TemperatureK, PressureBar),
                                density_superfluid(TemperatureK, PressureBar),
                                meanfreepath_adjustment_superfluid(TemperatureK, PressureBar), ballistic_switch)

    return res_width


def temperature_superfluid_from_width_obj(ObjOscillator, currentWidthHz, PressureBar):
    """
-    Calculates the temperature of the cell based on the current width for a given oscillator.
-
-    Made from old fortran code by Bradley, parameters are tuned for the Vibrating wires
    Much more sophistificated than a standard Stokes model
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
-    TemperatureK: temperature of helium mixture at SVP
-    """

    # define required temperature range, the module will not work correctly above ~100mK
    TemperatureK = np.linspace(0.1*temperature_critical_superfluid(PressureBar),temperature_critical_superfluid(PressureBar),num=1000) # generate temperature array
    # Calculate resonance width for all temperatures
    res_width = resonance_width_superfluid_obj(ObjOscillator, TemperatureK, PressureBar)

    # change float into numpy array for function len() to work
    try:
        len(currentWidthHz)
    except:
        currentWidthHz = np.array([currentWidthHz])
    currentWidthHz = np.array(currentWidthHz)

    currentTemperatureK = np.tile(np.NaN,len(currentWidthHz)) # array of NaN

    for tempindex in range(len(currentTemperatureK)):
        # Find the closest width an hence the temperature
        minindex = np.argmin(np.abs(res_width - np.tile(currentWidthHz[tempindex], len(TemperatureK))))
        currentTemperatureK[tempindex] = TemperatureK[minindex]

    return currentTemperatureK


def width_superfluid_ballistic_from_temperature(TemperatureK, coeff_He3B_Twidth, gapcoefficient = gap_coeff,
                                                PressureBar = 0.0):
    """
-    Calculates the thermalwidth an oscillator in the ballistic regime of He3 cell
-
     TemperatureK: temperature of helium-3 in K
     critical_temperature: critical temperature of superfluid transition [K]
     coeff_He3B_Twidth: coefficent for the wire (fork), size and material dependent
     gapcoeffcient: superfluid gap 1.76 according to BCS theory
-    return: thermal width of the oscillator
-    """

    return coeff_He3B_Twidth * np.exp(-gapcoefficient * temperature_critical_superfluid(PressureBar) / TemperatureK )

def temperature_superfluid_ballistic_from_width(thermalwidth, coeff_He3B_Twidth, gapcoefficient = gap_coeff,
                                                PressureBar = 0.0):
    """
-    Calculates the temperature of He3 cell in the ballistic regime only based on the current width for an oscillator.
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
     thermalwidth: thermal width of the oscillator
     critical_temperature: critical temperature of superfluid transition [K]
     coeff_He3B_Twidth: coefficent for the wire (fork), size and material dependent
     gapcoeffcient: superfluid gap 1.76 according to BCS theory
-    return: temperature of helium-3 in K
-    """

    return (-gapcoefficient * temperature_critical_superfluid(PressureBar)) / np.log(thermalwidth / coeff_He3B_Twidth)




def temperature_superfluid_ballistic_widthcoeff(wirewidth, wirecrossection, wiredensity,
                                                PressureBar = 0.0, gamma_coefficient = gamma_coefficient):
    """
-    Calculates the width coefficient for the vibrating object in the ballistic regime
-
-    wirewidth: wire or tuning fork width
     wirecrossection: crossection of the oscillator
     wiredensity: density of the oscillator
     PressureBar: Pressure in the cell
     gamma_coefficient: 0.28 for vibrating wires with 4.5 and 13.5 micron diameter
-    return: temperature width coefficient
-    """

    return (gamma_coefficient * 8 * Fermi_momentum(PressureBar)**4 / Plank_const**3 * wirewidth / (wirecrossection * wiredensity))

def temperature_superfluid_ballistic_widthcoeff_obj(objOscillator, PressureBar = 0.0, gamma_coefficient = gamma_coefficient):
    """
-    Calculates the temperature of He3 cell in the ballistic regime only based on the current width for an oscillator.
-
-    ObjOscillator: Oscillator object (wire or tuning fork) from mod_oscillator
     PressureBar: Pressure in the cell
     gamma_coefficient: 0.28 for vibrating wires with 4.5 and 13.5 micron diameter
-    return: temperature width coefficient
-    """
    if objOscillator.coeff_gamma == None:
        temp_gamma = gamma_coefficient
    else:
        temp_gamma = objOscillator.coeff_gamma

    return temperature_superfluid_ballistic_widthcoeff(objOscillator.width(), objOscillator.crossectionarea(),
                        objOscillator.density(), PressureBar = PressureBar, gamma_coefficient = temp_gamma)

def reduction_coeffcient(PressureBar = 0.0):
    """
-    Calculates the reduction coeffcient (Fermi momentum / Boltzmann constant)
-
     PressureBar: Pressure in the cell
-    return: reduction coefficient
-    """

    return Fermi_momentum(PressureBar=PressureBar) / Boltzmann_const

def reduction_factor(TemperatureK, PressureBar = 0.0):
    """
-    Calculates the reduction factor (Fermi momentum / (Boltzmann_constant*TemperatureK)
-
     PressureBar: Pressure in the cell
-    return: reduction factor
-    """

    return reduction_coeffcient(PressureBar=PressureBar) / TemperatureK

def width_thermal_ballistic_nonlinear_correction(widththermal, velocity, oscillator_lambda, TemperatureK, PressureBar = 0.0):
    """
-    Calculates the correct (low velocity) thermal width that appears smaller due to the Andreev reflection of quasiparticles
-
     widththermal: measured thermal width of the wire at a given velocity
     velocity: given velocity
     oscillator_lambda: for the wire approximately 0.8
     TemperatureK: cell temperature
     PressureBar: Pressure in the cell
-    return: corrected width
-    """

    nonlinear_correction = oscillator_lambda * reduction_factor(TemperatureK, PressureBar=PressureBar) * velocity / \
          (1.0 - np.exp(-oscillator_lambda * reduction_factor(TemperatureK, PressureBar=PressureBar) * velocity))
    return widththermal * nonlinear_correction

def width_thermal_ballistic_nonlinear_correction_obj(objObject, widththermal, velocity, gapcoefficient = gap_coeff, PressureBar = 0.0, iterations = 3):
    """
-    Calculates the correct (low velocity) thermal width from the smaller width due to the Andreev reflection of quasiparticles
-
     objObject: class corresponding to a given oscillator object
     widththermal: measured thermal width of the wire at a given velocity
     velocity: given velocity
     gapcoefficient: value of superfluid gap 1.76Tc by default
     PressureBar: Pressure in the cell
-    return: corrected width
-    """

    tempwidth = widththermal
    for tempintex in range(iterations):
        TemperatureK = temperature_superfluid_ballistic_from_width(tempwidth, objObject.coeff_He3B_Twidth,
                                    gapcoefficient = gapcoefficient, PressureBar = PressureBar)

        tempwidth = width_thermal_ballistic_nonlinear_correction(widththermal, velocity,
                                    objObject.coeff_lambda, TemperatureK, PressureBar = 0.0)

    return tempwidth


def width_thermal_ballistic_nonlinear_correction_reversed(widththermal, velocity, oscillator_lambda, TemperatureK, PressureBar = 0.0):
    """
-    Calculates the expected thermal width at high velocity from the low velocity (true) width due to the Andreev reflection of quasiparticles
-    needs correct temperature to be given
-
     widththermal: measured thermal width of the wire at a given velocity
     velocity: given velocity
     oscillator_lambda: for the wire approximately 0.8
     TemperatureK: cell temperature
     PressureBar: Pressure in the cell
-    return: corrected width
-    """

    nonlinear_correction = oscillator_lambda * reduction_factor(TemperatureK, PressureBar=PressureBar) * velocity / \
          (1.0 - np.exp(-oscillator_lambda * reduction_factor(TemperatureK, PressureBar=PressureBar) * velocity))
    return widththermal / nonlinear_correction


def widthparameter_ballistic(widththermal, TemperatureK, gapcoefficient = gap_coeff, PressureBar = 0.0):
    """
-    Calculates the width parameter: Width parameter = T[K]*Width[Hz]*( 1.76Tc + T )
-
     widththermal: measured thermal width of the wire
     TemperatureK: cell temperature
     PressureBar: Pressure in the cell
-    return: width_parameter
-    """

    widthparameter = TemperatureK * widththermal * (gapcoefficient * temperature_critical_superfluid(PressureBar) + TemperatureK)
    return widthparameter

def widthparameter_ballistic_obj(objObject, widththermal, gapcoefficient = gap_coeff, PressureBar = 0.0):
    """
-    Calculates the width parameter for an object from the thermal width: Width parameter = T[K]*Width[Hz]*( 1.76Tc + T )
     Determines the temperature from a given thermal width
-
     objObject: is a wire, tuning fork or other oscillator
     widththermal: measured thermal width of the wire
     PressureBar: Pressure in the cell
-    return: width_parameter
-    """
    TemperatureK = temperature_superfluid_ballistic_from_width(widththermal, objObject.coeff_He3B_Twidth, gapcoefficient = gapcoefficient, PressureBar = PressureBar)
    widthparameter = TemperatureK * widththermal * (gapcoefficient * temperature_critical_superfluid(PressureBar) + TemperatureK)
    return widthparameter

def width_thermal_effective_ballistic_obj(objObject, widththermal, objReference):
    """
-    Calculates the effective thermal width for an object based on referenced resonator
     Example: m1 (13.5 micron wire) compared with mmm3 (4.5 micron wire)
     the real thermal width of m1 needs to be divided by ~3 to be compared with thermal width of mmm3
-
     objObject: is a wire, tuning fork or other oscillator
     widththermal: measured thermal width of the wire
     objReference: is a wire, tuning fork or other oscillator to whos sensitivity width is recalculated
-    return: width_thermal_effective
-    """

    width_thermal_effective = widththermal * objReference.coeff_He3B_Twidth / objObject.coeff_He3B_Twidth
    return width_thermal_effective


############## DAQ table functions #####################

def daq_width_nonlinear_correction_ballistic_bulkT(tabDAQ, list_devices, list_channels, str_widthtrue = 'widthtrue',
        str_temperaturebulk ='TW29temperature', str_velocity = 'velocity', str_width='widththerm', slice_beg = None,
        slice_end = None, str_channel_header = 'TW', str_daqonly_header = 'DAC', PressureBar = 0.0,
        printdebug = False):

    if printdebug:
        print('Calculating nonlinear correction to the resonators width based on bulk temperature. DAQ contains columns x rows: %i x %i'
              %(len(tabDAQ.columns), len(tabDAQ)))

    if len(list_channels) == 0:
        #determine channels measured with lock-ins
        nondaqchannels = fio.daq_measured_channels_nondaq(tabDAQ, str_daqonly_header = str_daqonly_header)
    else:
        nondaqchannels = list_channels

    for tempindex in range(len(nondaqchannels)):
        tabDAQ.loc[slice_beg:slice_end,"%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_widthtrue)] = \
            width_thermal_ballistic_nonlinear_correction( #widththermal, velocity, oscillator_lambda, TemperatureK, PressureBar = 0.0
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_width)],
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_velocity)],
                list_devices[nondaqchannels[tempindex]-1].oscillator.coeff_lambda,
                tabDAQ.loc[slice_beg:slice_end, str_temperaturebulk],
                PressureBar = PressureBar)

    if printdebug:
        if len(nondaqchannels):
            print('current nonDAQ channel columns: "%s"'%(fio.daq_channel_column_names(tabDAQ, nondaqchannels[-1],
                                                                        str_channel_header = str_channel_header)))
        print(' ')

    return(tabDAQ)


def daq_width_nonlinear_correction_ballistic_obj(tabDAQ, list_devices, list_channels, str_widthtrue = 'widthtrue',
        str_velocity = 'velocity', str_width='widththerm', slice_beg = None, slice_end = None, str_channel_header = 'TW',
        str_daqonly_header = 'DAC', PressureBar = 0.0, iterations = 3, printdebug = False):

    if printdebug:
        print('Calculating true resonator width from velocity and lambda. DAQ contains columns x rows: %i x %i'
        %(len(tabDAQ.columns), len(tabDAQ)))

    if len(list_channels) == 0:
        #determine channels measured with lock-ins
        nondaqchannels = fio.daq_measured_channels_nondaq(tabDAQ, str_daqonly_header = str_daqonly_header)
    else:
        nondaqchannels = list_channels

    for tempindex in range(len(nondaqchannels)):
        tabDAQ.loc[slice_beg:slice_end,"%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_widthtrue)] = \
            width_thermal_ballistic_nonlinear_correction_obj(#objObject, widththermal, velocity, gapcoefficient = gap_coeff, PressureBar = 0.0, iterations = 3)
                list_devices[nondaqchannels[tempindex]-1].oscillator,
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_width)],
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_velocity)],
                PressureBar = PressureBar, iterations=iterations)

    if printdebug:
        if len(nondaqchannels):
            print('current nonDAQ channel columns: "%s"'%(fio.daq_channel_column_names(tabDAQ, nondaqchannels[-1],
                                                                        str_channel_header = str_channel_header)))
        print(' ')

    return(tabDAQ)



def daq_width_nonlinear_correction_ballistic(tabDAQ, list_devices, list_channels, str_widthtrue = 'widthtrue',
        str_temperature ='temperature', str_velocity = 'velocity', str_width='widththerm', slice_beg = None,
        slice_end = None, str_channel_header = 'TW', str_daqonly_header = 'DAC', PressureBar = 0.0,
        printdebug = False):

    if printdebug:
        print('Calculating nonlinear correction to the resonators width. DAQ contains columns x rows: %i x %i'
        %(len(tabDAQ.columns), len(tabDAQ)))

    if len(list_channels) == 0:
        #determine channels measured with lock-ins
        nondaqchannels = fio.daq_measured_channels_nondaq(tabDAQ, str_daqonly_header = str_daqonly_header)
    else:
        nondaqchannels = list_channels

    for tempindex in range(len(nondaqchannels)):
        tabDAQ.loc[slice_beg:slice_end,"%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_widthtrue)] = \
            width_thermal_ballistic_nonlinear_correction( #widththermal, velocity, oscillator_lambda, TemperatureK, PressureBar = 0.0
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_width)],
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_velocity)],
                list_devices[nondaqchannels[tempindex]-1].oscillator.coeff_lambda,
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_temperature)],
                PressureBar = PressureBar)

    if printdebug:
        if len(nondaqchannels):
            print('current nonDAQ channel columns: "%s"'%(fio.daq_channel_column_names(tabDAQ, nondaqchannels[-1],
                                                                        str_channel_header = str_channel_header)))
        print(' ')

    return(tabDAQ)


def daq_temperature_superfluid_ballistic(tabDAQ, list_devices, list_channels = [], str_temperature ='temperature',
                                         str_width='widththerm', slice_beg = None, slice_end = None,
                                         str_channel_header = 'TW', str_daqonly_header = 'DAC',
                                         gapcoefficient = gap_coeff, PressureBar = 0.0):
# Calculates temperature of the wires for the DAQ
#
    if len(list_channels) == 0:
        #determine channels measured with lock-ins
        nondaqchannels = fio.daq_measured_channels_nondaq(tabDAQ, str_daqonly_header = str_daqonly_header)
    else:
        nondaqchannels = list_channels

    for tempindex in range(len(nondaqchannels)):
        tabDAQ.loc[slice_beg:slice_end,"%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_temperature)] = \
            temperature_superfluid_ballistic_from_width(
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_width)],
                list_devices[nondaqchannels[tempindex]-1].oscillator.coeff_He3B_Twidth, gapcoefficient = gapcoefficient, PressureBar = PressureBar)
    return(tabDAQ)

def daq_width_thermal_effective_ballistic(tabDAQ, list_devices, objReference, list_channels = [],
                                          str_width_effective ='widththerm_effective',
                                          str_width='widththerm', slice_beg = None, slice_end = None,
                                          str_channel_header = 'TW', str_daqonly_header = 'DAC',
                                          printdebug = False):

# Calculates the effective thermal width of the channels (devices = object + surrounding) wires for the DAQ
# uses reference object, hance without surrounding. To use the device as reference, for example ch_mmm1, set to ch_mmm1.oscillator
#
    if len(list_channels) == 0:
        #determine channels measured with lock-ins
        nondaqchannels = fio.daq_measured_channels_nondaq(tabDAQ, str_daqonly_header = str_daqonly_header)
    else:
        nondaqchannels = list_channels

    for tempindex in range(len(nondaqchannels)):
        tabDAQ.loc[slice_beg:slice_end,"%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_width_effective)] = \
            width_thermal_effective_ballistic_obj(
                list_devices[nondaqchannels[tempindex]-1].oscillator,
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_width)],
                objReference)

    if printdebug:
        if len(nondaqchannels):
            print('current nonDAQ channel columns: "%s"'%(fio.daq_channel_column_names(tabDAQ, nondaqchannels[-1],
                                                                        str_channel_header = str_channel_header)))
        print(' ')

    return(tabDAQ)



def daq_widthparameter_ballistic(tabDAQ, list_devices=[], list_channels = [], str_widthpar ='widthpar', str_temperature ='temperature',
                                         str_width='widththerm', slice_beg = None, slice_end = None,
                                         str_channel_header = 'TW', str_daqonly_header = 'DAC',
                                         gapcoefficient = gap_coeff, PressureBar = 0.0):

    if len(list_channels) == 0:
        #determine channels measured with lock-ins
        nondaqchannels = fio.daq_measured_channels_nondaq(tabDAQ, str_daqonly_header = str_daqonly_header)
    else:
        nondaqchannels = list_channels

    for tempindex in range(len(nondaqchannels)):
        tabDAQ.loc[slice_beg:slice_end,"%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_widthpar)] = \
            widthparameter_ballistic(
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_width)],
                tabDAQ.loc[slice_beg:slice_end, "%s%i%s" % (str_channel_header, nondaqchannels[tempindex], str_temperature)],
                gapcoefficient = gapcoefficient, PressureBar = PressureBar)
    return(tabDAQ)

#####################################################################
####  Functions to work with ballistic beams of quasiparticles ######
#####################################################################

def quasiparticle_xip(quasiparticle_momentum, pressure = 0.0):

    quasiparticle_xip = quasiparticle_momentum[0]**2 + quasiparticle_momentum[1]**2 + quasiparticle_momentum[2]**2
    return(quasiparticle_xip)

def sph2cart(azimuth_theta, elevation_phi, radius_r, phi_zero_on_xyplane = True, round_2_zero = None, printdebug = False):
#SPH2CART Transform spherical to Cartesian coordinates.
#   [X,Y,Z] = SPH2CART(TH,PHI,R) transforms corresponding elements of
#   data stored in spherical coordinates (azimuth TH, elevation PHI,
#   radius R) to Cartesian coordinates X,Y,Z.  The arrays TH, PHI, and
#   R must be the same size (or any of them can be scalar).  TH and
#   PHI must be in radians.
#
#   TH is the counterclockwise angle in the xy plane measured from the
#   positive x axis.
#   PHI is the elevation angle from the xy plane.
#   round_2_zero: None or number below which rounding to zero takes place

    if phi_zero_on_xyplane:
        elevation_phi = elevation_phi
    else:
        elevation_phi = 0.5*np.pi - elevation_phi

    cartesian_z = radius_r * np.sin(elevation_phi)
    rcoselev = radius_r * np.cos(elevation_phi)
    cartesian_x = rcoselev * np.cos(azimuth_theta)
    cartesian_y = rcoselev * np.sin(azimuth_theta)

    if round_2_zero is not None:
        cartesian_x = np.where(np.less(np.absolute(cartesian_x), round_2_zero), +0.0e0, cartesian_x)
        cartesian_y = np.where(np.less(np.absolute(cartesian_y), round_2_zero), +0.0e0, cartesian_y)
        cartesian_z = np.where(np.less(np.absolute(cartesian_z), round_2_zero), +0.0e0, cartesian_z)

    return (cartesian_x, cartesian_y, cartesian_z)



def cart2sph(cartesian_x, cartesian_y, cartesian_z, phi_zero_on_xyplane = True, remove_NaN = False, alter_zero = 1.0e-38,
             round_2_zero = None):
#CART2SPH Transform  Cartesian coordinates to spherical.
#   [TH,PHI,R] = CART2SPH(X,Y,Z) transforms corresponding elements of
#   data stored in Cartesian coordinates X,Y,Z to spherical coordinates (azimuth TH, elevation PHI,
# #   radius R).  The arrays TH, PHI, and
#   R must be the same size (or any of them can be scalar).  TH and
#   PHI must be in radians.
#
#   TH is the counterclockwise angle in the xy plane measured from the
#   positive x axis.
#   PHI is the elevation angle from the xy plane.
#   round_2_zero: None or number below which rounding to zero takes place

    radius_r = np.sqrt(cartesian_x**2 + cartesian_y**2 + cartesian_z**2)
    if remove_NaN:
        cartesian_z = np.where(cartesian_z==0.0, cartesian_z + alter_zero, cartesian_z)
        cartesian_x = np.where(cartesian_x==0.0, cartesian_x + alter_zero, cartesian_x)
        cartesian_y = np.where(cartesian_y==0.0, cartesian_y + alter_zero, cartesian_y)

    azimuth_theta = np.arctan2(cartesian_y, cartesian_x)

    if phi_zero_on_xyplane:
        elevation_phi = np.arctan2(cartesian_z, np.sqrt(cartesian_x**2 + cartesian_y**2))
    else:
        elevation_phi = np.arctan2(np.sqrt(cartesian_x**2 + cartesian_y**2), cartesian_z)

    if round_2_zero is not None:
        azimuth_theta = np.where(np.less(np.absolute(azimuth_theta), round_2_zero), 0.0, azimuth_theta)
        elevation_phi = np.where(np.less(np.absolute(elevation_phi), round_2_zero), 0.0, elevation_phi)
        radius_r = np.where(np.less(np.absolute(radius_r), round_2_zero), 0.0, radius_r)

    return (azimuth_theta, elevation_phi, radius_r)


def polar2cart(azimuth_theta, radius_r):
#
#   azimuth_theta[rad] is the counterclockwise angle in the xy plane measured from the
#   positive x axis.
#
#   Class support for inputs TH,PHI,R:
#      float: double, single

    cartesian_x = radius_r * np.cos(azimuth_theta)
    cartesian_y = radius_r * np.sin(azimuth_theta)

    return (cartesian_x, cartesian_y)


def cartesian2cartesian_rotate90degrees_around_xaxis(x_coord, y_coord, z_coord, round_2_zero = None):
#   rotate cartesian coordinates around x-axis (x -> x, y->-z, z->y)
#   used for transforming particle on the orifice wall to normal surface of orifice
#   round_2_zero: None or number below which rounding to zero takes place

    x_coord_prime = x_coord
    y_coord_prime = -1.0 * z_coord
    z_coord_prime = y_coord

    if round_2_zero is not None:
        x_coord_prime = np.where(np.less(np.absolute(x_coord_prime), round_2_zero), +0.0e0, x_coord_prime)
        y_coord_prime = np.where(np.less(np.absolute(y_coord_prime), round_2_zero), +0.0e0, y_coord_prime)
        z_coord_prime = np.where(np.less(np.absolute(z_coord_prime), round_2_zero), +0.0e0, z_coord_prime)

    return (x_coord_prime, y_coord_prime, z_coord_prime)


"""
def spherical2spherical_rotate90degrees_around_yaxis(azimuth_theta, elevation_phi):
#
#   rotate spherical coordinates around y-axis (x->z, y->y, z->-x)
#   used for transforming particle on the orifice wall to normal surface of orifice
#   azimuth_theta[rad] is the counterclockwise angle in the xy plane measured from the
#   positive x axis.

    azimuth_theta_prime = np.arccos(-np.sin(elevation_phi)*np.cos(azimuth_theta))
    elevation_phi_prime = np.arctan( np.sin(elevation_phi)*np.tan(azimuth_theta))

    return (azimuth_theta_prime, elevation_phi_prime)

def spherical2spherical_rotate90degrees_around_xaxis(azimuth_theta, elevation_phi):
#
#   rotate spherical coordinates around x-axis (x -> x, y->z, z->-y)
#   used for transforming particle on the orifice wall to normal surface of orifice
#   azimuth_theta[rad] is the counterclockwise angle in the xy plane measured from the
#   positive x axis.

    azimuth_theta_prime = np.arctan(np.cos(elevation_phi)/(np.sin(elevation_phi)*np.cos(azimuth_theta)))
    elevation_phi_prime = np.arccos(-np.sin(elevation_phi)*np.sin(azimuth_theta))

    return (azimuth_theta_prime, elevation_phi_prime)
    
def spherical2spherical_rotate90degrees_around_origin(azimuth_theta, elevation_phi):
#
#   rotate spherical coordinates around x-axis (x -> y, y->z, z->-x)
#   used for transforming particle on the orifice wall to normal surface of orifice
#   azimuth_theta[rad] is the counterclockwise angle in the xy plane measured from the
#   positive x axis.

    azimuth_theta_prime = np.arctan(-np.tan(elevation_phi)*np.cos(azimuth_theta))
    elevation_phi_prime = np.arccos(-np.sin(elevation_phi)*np.sin(azimuth_theta))

    return (azimuth_theta_prime, elevation_phi_prime)
"""

def spherical2spherical_rotate90degrees_around_xaxis(azimuth_theta, elevation_phi, phi_zero_on_xyplane = True,
                                                     remove_NaN = True, alter_zero = 1.0e-38, round_2_zero = 1.e-8):
#   rotate spherical coordinates around x-axis (x -> x, y->-z, z->y)
#   used for transforming particle on the orifice wall to normal surface of orifice
#   azimuth_theta[rad] is the counterclockwise angle in the xy plane measured from the
#   positive x axis,
#   elevation_phi is zero at xy-axis if phi_zero_on_xyplane == True otherwise zero on z-axis.
#   round_2_zero: None or number below which rounding to zero takes place


    x_coord, y_coord, z_coord = sph2cart(azimuth_theta, elevation_phi, radius_r = 1.0,
                                         phi_zero_on_xyplane = phi_zero_on_xyplane, round_2_zero = round_2_zero)

    x_coord_prime, y_coord_prime, z_coord_prime = cartesian2cartesian_rotate90degrees_around_xaxis(x_coord, y_coord, z_coord,
                                                                        round_2_zero = round_2_zero)

    azimuth_theta_prime, elevation_phi_prime, radius_r = cart2sph(x_coord_prime, y_coord_prime, z_coord_prime,
                                                                  phi_zero_on_xyplane = phi_zero_on_xyplane,
                                                                  remove_NaN = remove_NaN, alter_zero = alter_zero,
                                                                  round_2_zero = round_2_zero)

    # make angle "theta" always positive
    azimuth_theta_prime[azimuth_theta_prime < 0.0] = azimuth_theta_prime[azimuth_theta_prime < 0.0] + 2 * np.pi

    return (azimuth_theta_prime, elevation_phi_prime)


def line3D_xy_at_z_eq2_0_from_2points(x_coord_init, y_coord_init, z_coord_init, x_coord, y_coord, z_coord, at_z_coordinate):
    # return coordinates x and y at z="orifice_escape_zcoordinate" for a line in 3D
    # line is given by 2 points (x_coord_init, y_coord_init, z_coord_init) and (x_coord, y_coord, z_coord)
    #follows from equation of line in 3D (x-x1)/l = (y-y1)/m = (z-z1)/n, where (l,m,n) is a vector going through a point
    # (x1,y1,z1). l = x2 - x1, m = y2 - y1, n = z2 - z1

    x_at_z_given = (at_z_coordinate - z_coord_init) / (z_coord - z_coord_init) *(x_coord - x_coord_init) + x_coord_init
    y_at_z_given = (at_z_coordinate - z_coord_init) / (z_coord - z_coord_init) *(y_coord - y_coord_init) + y_coord_init

    return x_at_z_given, y_at_z_given


def line_value_at_point(x_coord, line_slope, line_intercept):
    # calculate the value of the line at a given point

    return line_slope * x_coord + line_intercept

def line_xvalue_at_point(y_coord, line_slope, line_intercept, remove_NaN = True, alter_slope = 1.0e-38):
    # calculate the value of the line at a given point
    # remove_NaN: if line_slope == 0.0 change the slope by adding alter_slope to line_slope

    if remove_NaN:
        line_slope = np.where(line_slope==0.0, line_slope + alter_slope, line_slope)

    return (y_coord - line_intercept) / line_slope



def line_slope_from2points(x_coord, y_coord, x_coord2 = 0.0, y_coord2 = 0.0, remove_NaN = True, alter_x = 1.0e-38):
    # calculate the slope based on two points
    # remove_NaN [if x_coord == x_coord2 change the second by adding alter_x to x_coord2]

    diff_x = x_coord2 - x_coord
    diff_y = y_coord2 - y_coord

    if remove_NaN:
        diff_x = np.where(diff_x==0.0, diff_x + alter_x, diff_x)

    return np.true_divide(diff_y, diff_x)

def line_polarangle_from2points(x_coord, y_coord, x_coord2 = 0.0, y_coord2 = 0.0, remove_NaN = True, alter_x = 1.0e-38):
    # calculate the angle of the line based on two points
    # remove_NaN [if x_coord == x_coord2 change the second by adding alter_x to x_coord2]
    # direction is from 2nd coordinate to the first

    diff_x = x_coord2 - x_coord
    diff_y = y_coord2 - y_coord

    if remove_NaN:
        diff_x = np.where(diff_x==0.0, diff_x + alter_x, diff_x)

    return np.arctan2(diff_y, diff_x)



def line_intercept_from2points(x_coord, y_coord, x_coord2 = 0.0, y_coord2 = 0.0, remove_NaN = True, alter_x = 1.0e-38):
    # calculate the intercept with origin based on two points

    line_slope = line_slope_from2points(x_coord, y_coord, x_coord2, y_coord2, remove_NaN = remove_NaN, alter_x = alter_x)
    return line_intercept_from_point_and_slope(x_coord, y_coord, line_slope)

def line_intercept_from_point_and_slope(x_coord, y_coord, line_slope):
    # calculate the intercept based on a point and slope of the line

    return (y_coord - x_coord * line_slope)



def interception_line_and_circle(line_slope, line_intercept = 0, radius_circle = 1, round_2_zero = None):
# calculates crossing points of a line and a circle
# returns both solutions
# not tested for non-crossing case

    temp_slope2 = line_slope ** 2
    x_coordinate = (-line_slope*line_intercept + np.sqrt(temp_slope2*line_intercept**2 - (temp_slope2 + 1)*(line_intercept**2 - radius_circle**2)))/\
                   (temp_slope2 + 1)

    y_coordinate = line_value_at_point(x_coordinate, line_slope, line_intercept)

    x2_coordinate = (-line_slope*line_intercept - np.sqrt(temp_slope2*line_intercept**2 - (temp_slope2 + 1)*(line_intercept**2 - radius_circle**2)))/\
                   (temp_slope2 + 1)

    y2_coordinate = line_value_at_point(x2_coordinate, line_slope, line_intercept)

    if round_2_zero is not None:
        x_coordinate = np.where(np.less(np.absolute(x_coordinate), round_2_zero), +0.0e0, x_coordinate)
        y_coordinate = np.where(np.less(np.absolute(y_coordinate), round_2_zero), +0.0e0, y_coordinate)
        x2_coordinate = np.where(np.less(np.absolute(x2_coordinate), round_2_zero), +0.0e0, x2_coordinate)
        y2_coordinate = np.where(np.less(np.absolute(y2_coordinate), round_2_zero), +0.0e0, y2_coordinate)


    return x_coordinate, y_coordinate, x2_coordinate, y2_coordinate


def interception_particle_and_circle(polar_angle, x_coord_init = 0.0, y_coord_init = 0.0, radius_circle = 1.0,
                                     round_2_zero = 1e-12, printdebug = False):
# calculates crossing points of a particle trajectory and a circle
# assumes particle is fired from inner part of cylinder,
# not tested for non-crossing case

    polar_angle = angle_within_2pi(polar_angle)

    line_slope = np.tan(polar_angle)
    line_intercept = line_intercept_from_point_and_slope(x_coord_init, y_coord_init, line_slope)
    if printdebug:
        print('particle2cirle: polar_angle: %.1f, line: slope %.1f, intercept: %.1f'%(polar_angle, line_slope, line_intercept))

    x_coordinate, y_coordinate, x2_coordinate, y2_coordinate = interception_line_and_circle(line_slope,
                                                        line_intercept = line_intercept, radius_circle = radius_circle, round_2_zero = round_2_zero)

    if printdebug:
        temp_angle = angle_within_2pi(np.arctan2(y_coordinate - y_coord_init, x_coordinate - x_coord_init))
        print('particle2cirle: 1st solution: x_coord %.1f, y_coord: %.1f, polar angle: %.1f'%(x_coordinate, y_coordinate,
                                                    180.0/np.pi*temp_angle))

    if printdebug:
        temp2_angle = angle_within_2pi(np.arctan2(y2_coordinate - y_coord_init, x2_coordinate - x_coord_init))
        print('particle2cirle: 2nd solution: x_coord %.1f, y_coord: %.1f, polar angle: %.1f'%(x2_coordinate, y2_coordinate,
                                                    180.0/np.pi*temp2_angle))

    # choose appropriate coordinate
    arr_ycoord_max = y_coordinate * np.greater_equal(y_coordinate, y2_coordinate) + y2_coordinate * np.greater_equal(y2_coordinate, y_coordinate)
    arr_ycoord_min = y_coordinate * np.greater_equal(y2_coordinate, y_coordinate) + y2_coordinate * np.greater_equal(y_coordinate, y2_coordinate)
    y_true = arr_ycoord_max * np.greater_equal(np.pi, polar_angle) + arr_ycoord_min * np.greater(polar_angle, np.pi)

    #calculate corresponding x coordinate from the line equation
    x_true = line_xvalue_at_point(y_true, line_slope, line_intercept)

    if printdebug:
        print('particle2cirle: x_true %.1f, y_true: %.1f'%(x_true, y_true))

    return x_true, y_true


def angle_between_2_particlepaths(polar_angle, polar_angle2, radius = 1.0, radius2 = 1.0):
# calculates angle between two directions
# returns only accute angle
# not tested for non-crossing case

    x1, y1 = polar2cart(polar_angle, radius)
    x2, y2 = polar2cart(polar_angle2, radius2)

    dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
    det = x1*y2 - y1*x2      # determinant
    angle = np.arctan2(det, dot)  # atan2(y, x) or atan2(sin, cos)

    return angle

def angle_particlepath_2_circle_normal(polar_angle, x_coord_init = 0.0, y_coord_init = 0.0, radius = 1.0, printdebug = False):
# calculates angle between particle path and normal of the cylinder
# returns only accute angle
# not tested for non-crossing case

    x_true, y_true = interception_particle_and_circle(polar_angle, x_coord_init = x_coord_init,
                                                      y_coord_init = y_coord_init, radius_circle = radius)
    if printdebug:
        print('angle2normal: x_init %.1f, y_init: %.1f, x_true %.1f, y_true: %.1f'%(x_coord_init, y_coord_init, x_true, y_true))

    #polar angle of normal to this point
    normal_angle = line_polarangle_from2points(0.0, 0.0, x_true, y_true)

    angle_between = angle_between_2_particlepaths(polar_angle, normal_angle, radius = radius, radius2 = radius)

    if printdebug:
        print('angle2normal: incoming angle %.1f, angle of the normal %.1f, in between angle: %.1f'%(180.0/np.pi*polar_angle,
                                                                                180.0/np.pi*normal_angle, 180.0/np.pi*angle_between))

    return normal_angle, angle_between, x_true, y_true


def particlepath_specular_reflection_oncircle(normal_angle, angle_between, printdebug = False):

    inverse_normal_angle = angle_within_2pi(np.pi+normal_angle)

    specular_angle = angle_to_normal_2_polar(inverse_normal_angle, angle_between, printdebug = printdebug)

    if printdebug:
        print('specular reflection: angle of the normal %.1f, in between angle: %.1f, inverse normal %.1f, reflected %.1f '%(
               180.0/np.pi*normal_angle, 180.0/np.pi*angle_between, 180.0/np.pi*inverse_normal_angle, 180.0/np.pi*specular_angle))

    return specular_angle


def angle_interception_between_2lines(line_slope, line_slope2, accute_only = True):
# calculates angle between two lines
# returns only accute angle
# not tested for non-crossing case

    if accute_only:
        tangle_angle = np.arctan((line_slope2 - line_slope)/(1 + line_slope2 * line_slope))
    else:
        tangle_angle = np.arctan2((line_slope2 - line_slope),(1 + line_slope2 * line_slope))

    return tangle_angle


def line_slope_reflected_about_another_line(mirror_line_slope = 1, to_mirror_line_slope = 0):
# calculates slope of a line
# returns only accute angle
# not tested for non-crossing case

    new_line_slope = (-to_mirror_line_slope + mirror_line_slope**2 * to_mirror_line_slope + 2* mirror_line_slope) /\
                                (1 + 2 * to_mirror_line_slope * mirror_line_slope - mirror_line_slope**2)

    return new_line_slope

def angle_change_by_pi(angle):
    angle = np.where(angle > np.pi, angle - np.pi, angle + np.pi)
    return angle


def angle_within_2pi(arr_angle):
    #angle theat has to be positive

    arr_angle = np.where(arr_angle < 0.0, arr_angle + 2*np.pi, arr_angle)
    arr_angle = np.where(arr_angle >= 2*np.pi, arr_angle - 2*np.pi, arr_angle)

    return arr_angle



def angle_to_normal_2_polar(normal_angle, extra_angle, printdebug = False):
    #calculate change to azimuth angle theta

    polar_angle = normal_angle + extra_angle

    #angle should be within 0 and 2pi
    polar_angle_modified = angle_within_2pi(polar_angle)
    if printdebug:
        print('angle2normal->polar: angle of the normal %.1f, angle change: %.1f, polar angle %.1f, polar 2pi: %.1f'
              %(180.0/np.pi*normal_angle, 180.0/np.pi*extra_angle, 180.0/np.pi*polar_angle, 180.0/np.pi*polar_angle_modified))

    return(polar_angle_modified)



def quasiparticle_flux_sph2cart_cameraplane(arr_theta, arr_phi, distance = 1.0, x_coord_init = 0.0, y_coord_init = 0.0,
                                            phi_zero_on_xyplane = True):
# calculates the coordinates of quasiparticles on a face of the camera plane for the flux of qusiparticles emitted towards camera
#   arr_theta[rad] - theta angle (0 to 2pi)
#   arr_phi[rad] - phi angle (0 to 1/2pi)
#   distance - distance to camera plane
#   x,y,z - coordinates on the place
# consistent with Matlab spherical coordinates and
# http://mathworld.wolfram.com/SpherePointPicking.html

    if phi_zero_on_xyplane:
        temp_radius = distance / np.tan(arr_phi)
    else:
        temp_radius = distance * np.tan(arr_phi)

    arr_x = temp_radius * np.cos(arr_theta) + x_coord_init
    arr_y = temp_radius * np.sin(arr_theta) + y_coord_init
    arr_z = np.tile(distance,len(arr_x))

    return arr_x, arr_y, arr_z

def quasiparticle_flux_select_pixel_particles(pixel_x, pixel_y, pixel_radius, arr_theta_source, arr_phi_source, distance_2_source = 1.0,
                                x_coord_init_source = 0.0, y_coord_init_source = 0.0, phi_zero_on_xyplane = True, printdebug = True,
                                str_pixel_name='pix', open_geometry = False, capture_distance = 1.0):
#Flux_select_pixel_particles returns a fraction of particles going through a pixel and corresponing angles
#of particles into half unit sphere
#   pixel_x - pixel location X coordinate
#   pixel_y - pixel location Y coordinate
#   pixel_radius - Radius of the camera pixel
#   Arr_theta - Array of theta angles
#   Arr_phi - Array of phi angles
#   distance - distance to camera
#   Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Arr_entered_theta - Array of theta angles that are captured by pixel
#   Arr_entered_phi - Array of phi angles that are captured by pixel
#   entered_fraction - fraction of points selected from initial quasiparticles
#   open_geometry = False - all entered particled captured / "True" back side is open
#   capture_distance = 1.0 - added to the distance to the source

    particles = len(arr_theta_source)

    #Find indices within the camera pixel
    entered_indices = particle_inside_circle(pixel_x, pixel_y, pixel_radius, arr_theta_source, arr_phi_source,
                        distance = distance_2_source, x_coord_init = x_coord_init_source, y_coord_init = y_coord_init_source,
                        phi_zero_on_xyplane = phi_zero_on_xyplane)

    entered_indices = np.array(entered_indices[0])
    entered_particles = len(entered_indices)

    # Output arrays within the camera pixel
    entered_fraction = np.float(entered_particles) / np.float(particles)

    if open_geometry:
        transmitted_indices = particle_inside_circle(pixel_x, pixel_y, pixel_radius, arr_theta_source[entered_indices],
                            arr_phi_source[entered_indices], distance = distance_2_source + capture_distance,
                            x_coord_init = x_coord_init_source[entered_indices], y_coord_init=y_coord_init_source[entered_indices],
                                                  phi_zero_on_xyplane = phi_zero_on_xyplane)

        transmitted_indices = np.array(transmitted_indices[0])

        captured_indices =  np.array(list(set(entered_indices)-set(entered_indices[transmitted_indices])))
        captured_fraction = np.float(len(captured_indices)) / np.float(particles)

        if printdebug:
            print('%s %s: %s: %d, %s: %d, %s: %d, %s: %.2f, %s: %.2f' % ('Pixel -', str_pixel_name, 'incoming', particles,
            'entered', entered_particles, 'captured', len(captured_indices),
            'fraction entered %', entered_fraction * 100,
            'fraction captured %', captured_fraction * 100))

        if len(captured_indices) != 0:
            arr_captured_theta = arr_theta_source[captured_indices]
            arr_captured_phi = arr_phi_source[captured_indices]
        else:
            arr_captured_theta = []
            arr_captured_phi = []

    else:
        arr_captured_theta = arr_theta_source[entered_indices] # entered and captured are identical
        arr_captured_phi = arr_phi_source[entered_indices] # entered and captured are identical
        captured_fraction = entered_fraction

        if printdebug:
            print('%s %s: %s: %d, %s: %d, %s: %.2f%%' % ('Pixel -', str_pixel_name, 'incoming', particles,
                                                       'entered/captured', entered_particles,
                                                       'fraction entered/captured %', entered_fraction * 100))


    return captured_fraction, arr_captured_theta, arr_captured_phi, entered_indices





def particle_distance_2_circle_origin(circle_x, circle_y, arr_theta, arr_phi, distance = 1.0,
                           x_coord_init = 0.0, y_coord_init = 0.0, phi_zero_on_xyplane = True):
    #Convert spherical coordinates into circle plane
    arr_x, arr_y, arr_z = quasiparticle_flux_sph2cart_cameraplane(arr_theta, arr_phi, distance = distance,
                    x_coord_init = x_coord_init, y_coord_init = y_coord_init, phi_zero_on_xyplane = phi_zero_on_xyplane)

    #Calculate distance to a circle origin
    distance_2_circle = np.sqrt((arr_x - circle_x)**2 + (arr_y - circle_y)**2)

    return distance_2_circle, arr_x, arr_y, arr_z


def particle_inside_circle(circle_x, circle_y, circle_radius, arr_theta, arr_phi, distance = 1.0,
                           x_coord_init = 0.0, y_coord_init = 0.0, phi_zero_on_xyplane = True):

    #Calculate distance to a circle origin
    distance2pixel, arr_x, arr_y, arr_z = particle_distance_2_circle_origin(circle_x, circle_y, arr_theta, arr_phi,
                                        distance = distance, x_coord_init = x_coord_init, y_coord_init = y_coord_init,
                                        phi_zero_on_xyplane = phi_zero_on_xyplane)

    #Find indices within the camera pixel
    temp_ind = np.nonzero(distance2pixel < circle_radius)

    return temp_ind


def particles_inside_halfcircle(circle_x, circle_y, circle_radius, arr_theta, arr_phi, distance = 1.0,
                                x_coord_init = 0.0, y_coord_init = 0.0, phi_zero_on_xyplane = True):

    #Calculate distance to a circle origin
    distance2pixel, arr_x, arr_y, arr_z = particle_distance_2_circle_origin(circle_x, circle_y, arr_theta, arr_phi,
                                        distance = distance, x_coord_init = x_coord_init, y_coord_init = y_coord_init,
                                        phi_zero_on_xyplane = phi_zero_on_xyplane)

    #Find indices within the circle and top part of it
    temp_ind = np.nonzero((distance2pixel < circle_radius) * (arr_y >= circle_y))

    return temp_ind

def particles_inside_wireloop(wireloop_centre_x, wireloop_centre_y, wireloop_radius, wireloop_thickness, arr_theta, arr_phi,
                              distance = 1.0, x_coord_init = 0.0, y_coord_init = 0.0, phi_zero_on_xyplane = True):

    #Calculate distance to a circle origin
    distance2pixel, arr_x, arr_y, arr_z = particle_distance_2_circle_origin(wireloop_centre_x, wireloop_centre_y, arr_theta,
                                        arr_phi, distance = distance, x_coord_init = x_coord_init, y_coord_init = y_coord_init,
                                        phi_zero_on_xyplane = phi_zero_on_xyplane)

    #Find indices within the circle and top part of it
    temp_ind = np.nonzero((distance2pixel <= wireloop_radius + wireloop_thickness) *
                          (distance2pixel >= wireloop_radius - wireloop_thickness) *
                          (arr_y >= wireloop_centre_y))

    return temp_ind



def quasiparticle_flux_select_wireloop_particles(wireloop_centre_x, wireloop_centre_y, wireloop_radius, wireloop_thickness,
                                                 arr_theta, arr_phi, distance = 1.0, printdebug=True, x_coord_init = 0.0,
                                                 y_coord_init = 0.0, phi_zero_on_xyplane = True):
#Flux_select_pixel_particles returns a fraction of particles going through a pixel and corresponing angles
#of particles into half unit sphere
#   pixel_x - pixel location X coordinate
#   pixel_y - pixel location Y coordinate
#   pixel_radius - Radius of the camera pixel
#   Arr_theta - Array of theta angles
#   Arr_phi - Array of phi angles
#   distance - distance to camera
#%%Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Arr_selected_theta - Array of theta angles that go thru pixel
#   Arr_selected_phi - Array of phi angles that go thru pixel
#   %selected_fraction - fraction of points selected from initial quasiparticles
# consistent with Matlab spherical coordinates and
# http://mathworld.wolfram.com/SpherePointPicking.html

#    temp_ind = particles_inside_halfcircle(wire_x, wire_x, wire_radius, arr_theta, arr_phi, distance = 1.0)
    selected_indices = particles_inside_wireloop(wireloop_centre_x, wireloop_centre_y, wireloop_radius, wireloop_thickness,
                                         arr_theta, arr_phi, distance = distance, x_coord_init = x_coord_init,
                                         y_coord_init = y_coord_init, phi_zero_on_xyplane = phi_zero_on_xyplane)

    #Output arrays within the camera pixel
    arr_selected_theta = arr_theta[selected_indices]
    arr_selected_phi = arr_phi[selected_indices]
    selected_fraction = np.float(len(arr_selected_theta))/np.float(len(arr_theta))

    if printdebug:
        print('%s: %d %s: %d %s: %.4f'%('given particles', len(arr_theta), 'selected particles on wire loop', len(arr_selected_theta),
                                    'fraction selected %', selected_fraction*100))

    return selected_fraction, arr_selected_theta, arr_selected_phi, selected_indices


def quasiparticle_flux_into_solidangle(particles=1e5, max_elevation_phi=np.pi/2, str_distribution='uniform', phi_zero_on_xyplane=False):
# Flux_GeneratedAngles_ParticlesUniformHalfSphere returns angles for uniform distibution ...
#of particles into half unit sphere
#   particles - number of particles
#   random_theta[rad] - theta angle (0 to 2pi)
#   random_phi[rad] - phi angle (0 to 1/2pi)
# consistent with Matlab spherical coordinates and
# http://mathworld.wolfram.com/SpherePointPicking.html

# needs updating with phi zero on the plane

    arr_random_azimuth_theta = 2 * np.pi * np.random.rand(particles)

    if phi_zero_on_xyplane:
        arr_random_elevation_phi = np.tile(np.NAN, particles)
    else:
        if str_distribution == 'uniform':
            arr_random_elevation_phi = np.arccos(np.random.rand(particles)*(1-np.cos(max_elevation_phi)) + np.cos(max_elevation_phi))
        elif str_distribution == 'cos_theta':
            arr_random_elevation_phi = np.arccos(np.sqrt(np.random.rand(particles))*(1-np.cos(max_elevation_phi)) + np.cos(max_elevation_phi))
        else:
            arr_random_elevation_phi = np.tile(np.NAN, particles)

    return (arr_random_azimuth_theta, arr_random_elevation_phi)


def quasiparticle_orifice_distribution(particles = 1e5, orifice_radius = 1, random_theta = 2 * np.pi,
                                            str_distribution = 'uniform'):
# quasiparticle_seed_distribution_orifice - generates quasiparticles on the orifice 2d positions

# in polar coordinates
# dA = \rho d\rho d\theta = 1/2 d\rho^2 d\theta
#   particles - number of particles
#   random_theta[rad] - theta angle (0 to 2pi) smaller angle than 2pi will produce a segment
#   orifice_radius - orifice radius

    arr_random_polar_theta = random_theta * np.random.rand(particles)

    if str_distribution == 'uniform':
        arr_random_radius = np.sqrt(np.random.rand(particles))*orifice_radius
    else:
        arr_random_radius = np.tile(np.NAN, particles)

    return (arr_random_polar_theta, arr_random_radius)



def quasiparticle_intersection_with_cylinder(arr_xcoord_initial, arr_ycoord_initial, arr_zcoord_initial,
                                        arr_azimuth_theta_initial, arr_elevation_phi_initial, printdebug = False):
# calculate 2D projection of the quasiparticle emitted from orifice or other plane in the cylinder
# towards the wall of cylinder polar angle and radius and arimuth angle

    #determine x, y coordinates on the cylinder wall
    normal_angle, angle_between, arr_xcoord, arr_ycoord = angle_particlepath_2_circle_normal(arr_azimuth_theta_initial,
                                                x_coord_init = arr_xcoord_initial, y_coord_init = arr_ycoord_initial,
                                                printdebug = printdebug)

    #determine z coordinate
    arr_distance_inplane = np.sqrt((arr_xcoord - arr_xcoord_initial)**2 + (arr_ycoord - arr_ycoord_initial)**2)
    arr_zcoord = arr_zcoord_initial + arr_distance_inplane / (np.tan(arr_elevation_phi_initial))

    return (arr_xcoord, arr_ycoord, arr_zcoord, normal_angle, angle_between)


def quasiparticle_reflection_on_cylinder_specular(angle_of_normal, angle_between_incoming_and_normal,
                                        arr_elevation_phi_initial, printdebug = False):
# calculate 2D projection of the quasiparticle emitted from orifice or other plane in the cylinder
# towards the wall of cylinder polar angle and radius and arimuth angle

    # set new azimuth (elevation unchanged) angles after a specular reflection
    arr_theta = particlepath_specular_reflection_oncircle(angle_of_normal, angle_between_incoming_and_normal, printdebug = printdebug)

    #return elevation angle phi unchanged (specular reflection)
    arr_phi = arr_elevation_phi_initial

    return (arr_theta, arr_phi)

def quasiparticle_reflection_on_cylinder_diffuse(angle_of_normal, angle_between_incoming_and_normal,
                                        arr_elevation_phi_initial, printdebug = False, phi_zero_on_xyplane = True):
# calculate 2D projection of the quasiparticle emitted from orifice or other plane in the cylinder
# towards the wall of cylinder polar angle and radius and arimuth angle

    arr_random_azimuth_theta, arr_random_elevation_phi = quasiparticle_flux_into_solidangle(particles=np.shape(arr_elevation_phi_initial)[0])

    #angle between cylinder normal and outgoing particle
    theta_angle_to_normal, arr_phi = spherical2spherical_rotate90degrees_around_xaxis(arr_random_azimuth_theta, arr_random_elevation_phi, phi_zero_on_xyplane = phi_zero_on_xyplane)

    if printdebug and np.greater(np.shape(arr_elevation_phi_initial)[0], 0):
        print('diffuse scattering theta: min %.1f, max: %.1f, mean:, %.1f, phi:  min %.1f, max: %.1f, mean:, %.1f'%(np.min(theta_angle_to_normal),
                np.max(theta_angle_to_normal), np.mean(theta_angle_to_normal), np.min(arr_phi), np.max(arr_phi), np.mean(arr_phi)))

#    arr_theta = angle_to_normal_2_polar(inverse_normal_angle, theta_angle_to_normal, printdebug = printdebug)

    # set new azimuth angles after a reflection with newly generated random angles
    arr_theta = particlepath_specular_reflection_oncircle(angle_of_normal, theta_angle_to_normal - np.pi, printdebug = printdebug)

    return (arr_theta, arr_phi)


def quasiparticle_reflection_on_cylinder_mixed(angle_of_normal, angle_between_incoming_and_normal,
                                        arr_elevation_phi_initial, degree_of_specularity=1.0, printdebug = False, phi_zero_on_xyplane = False):
# calculate 2D projection of the quasiparticle emitted from orifice or other plane in the cylinder
# towards the wall of cylinder polar angle and radius and arimuth angle

# degree_of_specularity = 0 (1) fully diffuse (specular)
    particles_local = np.shape(arr_elevation_phi_initial)[0]
    #determine a maximum spread of the beam emitted after collision with the wall
    max_elevation_phi_local = np.pi/2 * (1.0 - degree_of_specularity)

    arr_random_azimuth_theta, arr_random_elevation_phi = quasiparticle_flux_into_solidangle(particles=particles_local,
                                    max_elevation_phi=max_elevation_phi_local, phi_zero_on_xyplane=phi_zero_on_xyplane)

    #angle between cylinder normal and outgoing particle
    theta_angle_to_normal, arr_phi = spherical2spherical_rotate90degrees_around_xaxis(arr_random_azimuth_theta,
                                                arr_random_elevation_phi, phi_zero_on_xyplane = phi_zero_on_xyplane)

    #work out the final angle for phi
    arr_phi_mixed = arr_phi * (1.0 - degree_of_specularity) + degree_of_specularity * arr_elevation_phi_initial

    if printdebug and np.greater(particles_local, 0):
        print('diffuse scattering theta: min %.1f, max: %.1f, mean:, %.1f, phi:  min %.1f, max: %.1f, mean:, %.1f'%(np.min(theta_angle_to_normal),
                np.max(theta_angle_to_normal), np.mean(theta_angle_to_normal), np.min(arr_phi_mixed), np.max(arr_phi_mixed), np.mean(arr_phi_mixed)))

    theta_angle_mixed = (theta_angle_to_normal) + degree_of_specularity * angle_between_incoming_and_normal

    # set new azimuth angles after a reflection with newly generated random angles
    arr_theta = particlepath_specular_reflection_oncircle(angle_of_normal, theta_angle_mixed, printdebug=printdebug)

    return (arr_theta, arr_phi_mixed)





def angle_qp_emission_cylinder_surface(velocity_wire, polar_angle = np.pi/2.0, PressureBar = 0.0, resulting_angle_in_degrees = False):
# BBR PRL Fisher based on alpha*p_F*v = \Delta - p_F * v * cos(phi)
    wire_radius = 1.0

    angle_arccos = Landau_velocity(PressureBar) / velocity_wire - osc.velocity_enhancement_potential_flow_around_cylinder(wire_radius, polar_angle, wire_radius)

#    if angle_arccos > 1.0: #for low velocities emission does not take place
#        angle_arccos = 1.0

    if angle_arccos < 0.0: #for high velocities emission angle does not exceed pi/2
        angle_arccos = 0.0

    angle = np.arccos(angle_arccos)
    if resulting_angle_in_degrees:
        angle = angle * 180.0 / np.pi

    return angle

def angle_qp_emission_from_cylinder_surface(velocity_wire, polar_angle = np.pi/2.0, PressureBar = 0.0, resulting_angle_in_degrees = False):
# take into account angles between p_f and velocity
    wire_radius = 1.0

    angle_arccos = Landau_velocity(PressureBar) / (velocity_wire * (np.sin(polar_angle) + osc.velocity_enhancement_potential_flow_around_cylinder(wire_radius, polar_angle, wire_radius)))

#    if angle_arccos > 1.0: #for low velocities emission does not take place
#        angle_arccos = 1.0

    if angle_arccos < 0.0: #for high velocities emission angle does not exceed pi/2
        angle_arccos = 0.0

    angle = np.arccos(angle_arccos)
    if resulting_angle_in_degrees:
        angle = angle * 180.0 / np.pi

    return angle



###########################################################################
#######   BBR functions ###########
###########################################################################


def time_constant_BBR_ballistic_superfluid(volume, orifice_area, group_velocity_quasiparticles):
# Approximate time costant of the black box radiator in superfluid He3
# volume - volume of the BBR box [m^3]
# orifice_area - hole area in the BBR [m^2]
# group_velocity_quasiparticles - [m/s] in ballistic regime approximately 1/4 of Fermi velocity

    time_costant = 4 * volume / (orifice_area * group_velocity_quasiparticles)

    return time_costant

def width_parameter_BBR_factor(orifice_area, thermometer_wire_diameter, thermometer_wire_mass_perunitlength,
                               gamma_coefficient = gamma_coefficient, PressureBar = 0.0, orifice_area_multiplier = 1.0,
                               units_kB = True):

    width_parameter_factor = gamma_coefficient * (2.0 * thermometer_wire_diameter * Fermi_momentum(PressureBar)**2) / \
                      (np.pi * thermometer_wire_mass_perunitlength * orifice_area_multiplier * orifice_area * Boltzmann_const)

    if units_kB:
        width_parameter_factor = width_parameter_factor / Boltzmann_const

    return width_parameter_factor



def width_parameter_BBR_calibration(power_appliedW, orifice_area = np.pi*0.15e-3**2, thermometer_wire_diameter = 4.5e-6,
                                    thermometer_wire_mass_perunitlength = 6.05e3*np.pi*4.5e-6**2/4.0,
                                    gamma_coefficient = gamma_coefficient, PressureBar = 0.0, orifice_area_multiplier = 1.0, units_kB = True):

    width_parameter = width_parameter_BBR_factor(orifice_area, thermometer_wire_diameter, thermometer_wire_mass_perunitlength,
                               gamma_coefficient = gamma_coefficient, PressureBar = PressureBar, orifice_area_multiplier = orifice_area_multiplier,
                               units_kB = units_kB) * power_appliedW

    if units_kB:
        width_parameter = width_parameter / Boltzmann_const

    return width_parameter

