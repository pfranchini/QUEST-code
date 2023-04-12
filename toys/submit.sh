#!/bin/bash

# Submit in nohup all the scripts for a full toy simulation (lockin, SQUID, shot-noise) given a config.py file

cp config.py output/
nohup python lockin_toy_energy.py > lockin_toy_energy.out &
nohup python lockin_toy.py > lockin_toy.out &
nohup python squid_toy_energy.py >squid_toy_energy.out &
nohup python squid_toy.py >squid_toy.out &
nohup python shot_toy_energy.py > shot_toy_energy.out &
nohup python shot_toy.py > shot_toy.out &

