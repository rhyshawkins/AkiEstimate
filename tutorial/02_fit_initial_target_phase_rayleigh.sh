#!/bin/bash

if [ -z "$1" ]
then
  STATIONPAIR=HOT05_HOT25
else
  STATIONPAIR=$1
fi


mkdir -p Initial_${STATIONPAIR}

FMIN=0.025
FMAX=0.35
RHO=0.5e3
VS=0.5e3
XI=0.05
VPVS=0.05

SKIP=15

#gdb --args 
/usr/bin/time -o Initial_${STATIONPAIR}/opt.time ../InitialPhase/optimizer/optimizerayleigh \
    -i ../example_data/RayleighResponse/dispersion_${STATIONPAIR}.txt \
    -C InitialPhase/phase_${STATIONPAIR}.rayleigh \
    -o Initial_${STATIONPAIR}/opt \
    -r ../Reference/models/greenetal+lidettrick_smooth_edgecorrected.txt \
    -f $FMIN \
    -F $FMAX \
    -R $RHO \
    -V $VS \
    -X $XI \
    -S $VPVS \
    -M 0 \
    -N 30 \
    -e 1.0 \
    -T $SKIP


