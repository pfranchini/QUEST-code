## Hybrid-Toy Simulation

### Description:
 - QUEST-DMC WP1 full toy simulation of a train of pulses from an energy pdf
 - Reads AM results from the WP1 Fridge 4 run present in `/data` to simulate a train of pulses given a specific run (`--time`) using the real noise from the run
 - Performs an analysis of the train of pulses to be compared with the data analysis, including:
    * noise
    * peak finding
    * resolution

### Usage:

#### Required input:
 - a config file (example in `config.py`)
 - an energy pdf: `energy|entries` (example in `energy_pdf.txt`)

```
python toy.py [-h] [--time TIME] [--config CONFIG] [--noise NOISE]

options:
  -h, --help       show this help message and exit
  --time TIME      Starting time [s]
  --config CONFIG  Config file
  --noise NOISE    Noise model (normal|shot|real)
```

#### E.g.
```
python toy.py --config config.py --noise real --time 1675908556
```

![image](https://github.com/user-attachments/assets/ebd6f64d-ec72-4f0f-93d6-9ed5c6537561)
