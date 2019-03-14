#!/bin/bash

STATIONPAIR=HOT05_HOT25
mkdir -p Initial_${STATIONPAIR}

FMIN=0.025
FMAX=0.35
RHO=0.5e3
VS=0.5e3
XI=0.05
VPVS=0.05

../InitialPhase/optimizer/optimizejoint \
    -i ../example_data/LoveResponse/dispersion_${STATIONPAIR}.txt \
    -I ../example_data/LoveResponse/dispersion_${STATIONPAIR}.txt \
    -c InitialPhase/phase_${STATIONPAIR}.love \
    -C InitialPhase/phase_${STATIONPAIR}.rayleigh \
    -o Initial_${STATIONPAIR}/opt \
    -r ../Reference/models/greenetal+lidettrick_smooth_edgecorrected.txt \
    -f $FMIN \
    -F $FMAX \
    -R $RHO \
    -V $VS \
    -X $XI \
    -S $VPVS \
    -M 0 -A 0.05 \
    -N 1 \
    -e 1.0


