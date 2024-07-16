## Hybrid-Toy Simulation

Usage:
```
toy.py [-h] [--time TIME] [--config CONFIG] [--noise NOISE]

options:
  -h, --help       show this help message and exit
  --time TIME      Starting time [s]
  --config CONFIG  Config file
  --noise NOISE    Noise model (normal|shot|real)
```

E.g.
```
python toy.py --config config.py --noise real --time 1675908556
```