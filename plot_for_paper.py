'''
Plots on the same plot the Error vs Pressure, Diameter, T/Tc for the
lock-in amplfier and the SQUID readout toys produce before and saved
in a text file.
'''

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'lines.linewidth': 3})

# Pressure:

with open('lockin_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:] #skip first line
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]


with open('squid_toy-error-pressure.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

plt.plot(x1,y1,label='Lock-in amplifier')
plt.plot(x2,y2,label='SQUID readout')
plt.xlabel('Pressure [bar]')
plt.ylabel('Error [%]')
#plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("Error_pressure-comparison.pdf")
plt.show()

# Diameter:

x1=np.array([])
x2=np.array([])
y1=np.array([])
y2=np.array([])

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

plt.plot(x1*1e9,y1,label='Lock-in amplifier')
plt.plot(x2*1e9,y2,label='SQUID readout')
plt.xlabel('Diameter [nm]')
plt.ylabel('Error [%]')
plt.yscale('log')
plt.legend()
plt.savefig("Error_diameter-comparison.pdf")
plt.show()

# T/Tc:

with open('lockin_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:] #skip first line
    x1 = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]


with open('squid_toy-error-temperature.txt') as f:
    lines = f.readlines()[1:]
    x2 = [float(line.split()[0]) for line in lines]
    y2 = [float(line.split()[1]) for line in lines]

plt.plot(x1,y1,label='Lock-in amplifier')
plt.plot(x2,y2,label='SQUID readout')
plt.xlabel('T/T$_c$')
plt.ylabel('Error [%]')
#plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("Error_temperature-comparison.pdf")
plt.show()


