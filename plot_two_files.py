import matplotlib.pyplot as plt

with open('lockin-error-5bar.txt') as f:
    lines = f.readlines()
    x1 = [line.split()[0] for line in lines]
    y1 = [line.split()[1] for line in lines]


with open('squid_toy-error.txt') as f:
    lines = f.readlines()
    x2 = [line.split()[0] for line in lines]
    y2 = [line.split()[1] for line in lines]


print(y1)
    
plt.plot(x1,y1,label='Lock-in amplifier')
plt.plot(x2,y2,label='SQUID readout')
plt.xlabel('Energy [eV]')
plt.ylabel('Error [%]')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("Error_comparison.pdf")
plt.show()

