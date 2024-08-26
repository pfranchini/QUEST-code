## Hybrid-Toy Simulation

### Description:
 - QUEST-DMC WP1 full toy simulation of a train of pulses from an energy pdf
 - Reads AM results from the WP1 Fridge 4 run in `/data` to simulate a train of pulses given a specific run (`--time`)
 - Analysis of the train of pulses:
    * noise
    * peak finding
    * resolution

### Usage:
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