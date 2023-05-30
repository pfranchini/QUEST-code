'''
Plots on the same plot the Error vs Pressure, Diameter, T/Tc for the
lock-in amplfier and the SQUID readout toys produce before and saved
in a text file.

Two plots:
  E | T/Tc
  d | pressure

Input: txt files produced by the toy scripts.
Usage:  
   cd output direcotry
   python ../../plot_for_paper_shot.py

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

### First plot:
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Energy:

with open('lockin-error.txt') as f:
    lines = f.readlines()[1:]
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]

with open('squid_toy-error.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

with open('shot-error.txt') as f:
    lines = f.readlines()[1:]
    x3 = [float(line.split()[0]) for line in lines]
    y3 = [float(line.split()[1]) for line in lines]

ax1.text(0.05, 0.90, 'a', transform=ax1.transAxes, fontsize=16, verticalalignment='top', horizontalalignment='left', weight='bold')
ax1.plot(x1,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
ax1.plot(x2,y2,label='SQUID readout',color='red',linestyle='--')
ax1.plot(x3,y3,label='QP shot noise',color='green',linestyle=':')
ax1.set_xlabel('Energy [eV]')
ax1.set_ylabel('Error [%]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()

# T/Tc:

with open('lockin_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:] #skip first line
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]

popt, pcov = curve_fit(heat_capacity,x1,y1)
A = popt

with open('squid_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

with open('shot_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:]
    x3 = [float(line.split()[0]) for line in lines]
    y3 = [float(line.split()[1]) for line in lines]

ax2.text(0.05, 0.90, 'b', transform=ax2.transAxes, fontsize=16, verticalalignment='top', horizontalalignment='left', weight='bold')
ax2.plot(x1,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
ax2.plot(x2,y2,label='SQUID readout',color='red',linestyle='--')
ax2.plot(x3,y3,label='QP shot noise',color='green',linestyle=':')    
ax2.plot(x1, heat_capacity(np.array(x1), A),linestyle=':', marker='',color="black")
ax2.set_xlabel('T/T$_c$')
ax2.set_ylabel('Error [%]')
#ax2.xscale('log')
ax2.set_yscale('log')
#ax2.legend()

plt.tight_layout()
plt.subplots_adjust(wspace=0.3)
plt.savefig("Error_energy-temperature_comparison.pdf")
plt.savefig("Error_energy-temperature_comparison.png")
plt.show()


### Second plot:
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Diameter:

x1=np.array([])
x2=np.array([])
x3=np.array([])
y1=np.array([])
y2=np.array([])
y3=np.array([])

with open('lockin_toy-error-diameter.txt') as f:
    lines = f.readlines()[1:] #skip first line
    for line in lines:
        x1 = np.append(x1, float(line.split()[0]))
        y1 = np.append(y1, float(line.split()[1]))

with open('squid_toy-error-diameter.txt') as f:
    lines = f.readlines()[1:]
    for line in lines:
        x2 = np.append(x2, float(line.split()[0]))
        y2 = np.append(y2, float(line.split()[1]))

with open('shot_toy-error-diameter.txt') as f:
    lines = f.readlines()[1:]
    for line in lines:
        x3 = np.append(x3, float(line.split()[0]))
        y3 = np.append(y3, float(line.split()[1]))

ax1.text(0.05, 0.90, 'a', transform=ax1.transAxes, fontsize=16, verticalalignment='top', horizontalalignment='left', weight='bold')
ax1.plot(x1*1e9,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
ax1.plot(x2*1e9,y2,label='SQUID readout',color='red',linestyle='--')
ax1.plot(x3*1e9,y3,label='QP shot noise',color='green',linestyle=':')
ax1.set_xlabel('Diameter [nm]')
ax1.set_ylabel('Error [%]')
ax1.set_yscale('log')
ax1.legend()

# Pressure:

with open('lockin_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:] #skip first line
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]

with open('squid_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

with open('shot_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:]
    x3 = [float(line.split()[0]) for line in lines]
    y3 = [float(line.split()[1]) for line in lines]

ax2.text(0.05, 0.90, 'b', transform=ax2.transAxes, fontsize=16, verticalalignment='top', horizontalalignment='left', weight='bold')
ax2.plot(x1,y1,label='Lock-in amplifier',color='dodgerblue',linestyle='-')
ax2.plot(x2,y2,label='SQUID readout',color='red',linestyle='--')
ax2.plot(x3,y3,label='QP shot noise',color='green',linestyle=':')
ax2.set_xlabel('Pressure [bar]')
ax2.set_ylabel('Error [%]')
#ax2.xscale('log')
ax2.set_yscale('log')
#ax2.legend()

plt.tight_layout()
plt.subplots_adjust(wspace=0.3)
plt.savefig("Error_diameter-pressure_comparison.pdf")
plt.savefig("Error_diameter-pressure_comparison.png")
plt.show()
