#!/bin/bash

if [ -z "$1" ]
then
  STATIONPAIR=HOT05_HOT25
else
  STATIONPAIR=$1
fi


mkdir -p InitialPhase

python2 ../InitialPhase/scripts/estimate_rayleigh_phase_amplitude.py -p ../example_data -s $STATIONPAIR -o InitialPhase/phase_$STATIONPAIR --noshow --filter 3
