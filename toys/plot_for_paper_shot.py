'''
Plots on the same plot the Error vs Pressure, Diameter, T/Tc for the
lock-in amplfier and the SQUID readout toys produce before and saved
in a text file.

Input: txt files produced by the toy scripts.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'lines.linewidth': 3})
plt.rcParams.update({'xtick.direction': 'in'})
plt.rcParams.update({'ytick.direction': 'in'})

#============================================================

def heat_capacity(T,A):
    #return A*g*g/T*np.exp(-g/T);
    return A*np.exp(-1.76/T);

#============================================================


# Energy:

with open('output/lockin-error.txt') as f:
    lines = f.readlines()[1:]
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]

with open('output/squid_toy-error.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

with open('output/shot-error.txt') as f:
    lines = f.readlines()[1:]
    x3 = [float(line.split()[0]) for line in lines]
    y3 = [float(line.split()[1]) for line in lines]

plt.plot(x1,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
plt.plot(x2,y2,label='SQUID readout',color='red',linestyle='--')
plt.plot(x3,y3,label='QP shot noise',color='green',linestyle=':')
plt.xlabel('Energy [eV]')
plt.ylabel('Error [%]')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("Error_comparison.pdf")
plt.savefig("Error_comparison.png")
plt.show()

# Pressure:

with open('output/lockin_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:] #skip first line
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]

with open('output/squid_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

with open('output/shot_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:]
    x3 = [float(line.split()[0]) for line in lines]
    y3 = [float(line.split()[1]) for line in lines]

plt.plot(x1,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
plt.plot(x2,y2,label='SQUID readout',color='red',linestyle='--')
plt.plot(x3,y3,label='QP shot noise',color='green',linestyle=':')
plt.xlabel('Pressure [bar]')
plt.ylabel('Error [%]')
#plt.xscale('log')
plt.yscale('log')
#plt.legend()
plt.savefig("Error_pressure-comparison.pdf")
plt.savefig("Error_pressure-comparison.png")
plt.show()

# Diameter:

x1=np.array([])
x2=np.array([])
x3=np.array([])
y1=np.array([])
y2=np.array([])
y3=np.array([])

with open('output/lockin_toy-error-diameter.txt') as f:
    lines = f.readlines()[1:] #skip first line
    for line in lines:
        x1 = np.append(x1, float(line.split()[0]))
        y1 = np.append(y1, float(line.split()[1]))

with open('output/squid_toy-error-diameter.txt') as f:
    lines = f.readlines()[1:]
    for line in lines:
        x2 = np.append(x2, float(line.split()[0]))
        y2 = np.append(y2, float(line.split()[1]))

with open('output/shot_toy-error-diameter.txt') as f:
    lines = f.readlines()[1:]
    for line in lines:
        x3 = np.append(x3, float(line.split()[0]))
        y3 = np.append(y3, float(line.split()[1]))

plt.plot(x1*1e9,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
plt.plot(x2*1e9,y2,label='SQUID readout',color='red',linestyle='--')
plt.plot(x3*1e9,y3,label='QP shot noise',color='green',linestyle=':')
plt.xlabel('Diameter [nm]')
plt.ylabel('Error [%]')
plt.yscale('log')
#plt.legend()
plt.savefig("Error_diameter-comparison.pdf")
plt.savefig("Error_diameter-comparison.png")
plt.show()

# T/Tc:

with open('output/lockin_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:] #skip first line
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]

popt, pcov = curve_fit(heat_capacity,x1,y1)
A = popt

with open('output/squid_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

with open('output/shot_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:]
    x3 = [float(line.split()[0]) for line in lines]
    y3 = [float(line.split()[1]) for line in lines]

plt.plot(x1,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
plt.plot(x2,y2,label='SQUID readout',color='red',linestyle='--')
plt.plot(x3,y3,label='QP shot noise',color='green',linestyle=':')    
#plt.plot(x1, heat_capacity(np.array(x1), A)*5e9,linestyle=':', marker='',color="black")
plt.xlabel('T/T$_c$')
plt.ylabel('Error [%]')
#plt.xscale('log')
plt.yscale('log')
#plt.legend()
plt.savefig("Error_temperature-comparison.pdf")
plt.savefig("Error_temperature-comparison.png")
plt.show()


