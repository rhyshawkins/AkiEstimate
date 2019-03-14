#!/bin/bash

mkdir -p reference

MODEL=../Reference/models/greenetal+lidettrick_smooth_edgecorrected.txt

FMIN=0.001
FMAX=0.52

../Reference/mkreferencelove -i $MODEL -S 8192 -H 2.0 -f $FMIN -F $FMAX -o reference/reference_love_fine.txt

../Reference/mkreferencerayleigh -i $MODEL -S 8192 -H 2.0 -f $FMIN -F $FMAX -o reference/reference_rayleigh_fine.txt
