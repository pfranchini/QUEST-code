
### Generates the plots present in the note:
```
[]$ python code.py --help

code.py -p <pressure [bar]> -d <diameter [m]>
```

### SQUID readout plots as in Lev's note:
```
[]% python3 squid.py
```

### Extract energy from Viktor's data:
```
[]$ python viktor.py
```

### Study of the pulse peak time dependency:
```
[]$ python pulse.py -p <pressure [bar]> -d <diameter [m]> -t <relative temperature> -b <decay constant t_b>
```

## Toys:

### Generates the error for the energy deposition measurement, for an energy range:
```
[]$ python lockin_toy_energy.py
[]$ squid_toy_energy.py
```

### Generates the error for diameter, pressure and temperature, for the fixed energy deposition measurement:
```
[]$ python lockin_toy.py
[]$ squid_toy.py
```


