#!/bin/bash

STATIONPAIR=HOT05_HOT25

mkdir -p InitialPhase

python2 ../InitialPhase/scripts/estimate_joint_phase_amplitude.py -p ../example_data -s $STATIONPAIR -o InitialPhase/phase_$STATIONPAIR --noshow --filter 3
