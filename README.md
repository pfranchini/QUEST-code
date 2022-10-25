
### Generates the plots present in the note:
```
[]$ python code.py --help

code.py -p <pressure [bar]> -d <diameter [m]>
```

### Generates the error for the energy deposition measurement, for an energy range:
```
[]$ python toy.py
```

SQUID readout:
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