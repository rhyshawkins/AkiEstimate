# AkiEstimate

Estimate dispersion from Ambient Noise Cross-Correlations using Aki's spectral formulation

## Introduction

This software estimates interstation dispersion for Love and Rayleigh surface waves using
Aki's spectral formulation. The method is detailed in the publication:

  Hawkins R. and Sambridge M., ``An adjoint technique for estimation of interstation phase and
  group dispersion from ambient noise cross-correlations``, BSSA, 2019

This software is scientific software and is intended to demonstrate the functionality outlined
in the manuscript reference above. 

## Authors

Rhys Hawkins

## License

This code is released under a GPL v3 License.

A portion of the public domain LAPACK codes obtained from netlib.org
are included in this repository as a convenience under the
forwardmodel/dggev subdirectory as is licensed by the authors under
a modified BSD license (see http://www.netlib.org/lapack/LICENSE.txt).

## Prequisites

### For compilation

* g++/gfortran or equivalent (versions 8.2.x were used during development)
* FFTW http://www.fftw.org/
* GSL https://www.gnu.org/software/gsl/

### For processing and plotting

* Python 2.7
* Matplotlib
* Numpy

## Compilation

You may need to edit the Makefiles in Reference, InitialPhase/optimizer
and Phase/optimizer to adjust compiler used.

The top level Makefile compiles all the necessary components so assuming
the prerequisites are available as listed above, compilation should be
a matter of issuing a `make' command from the top level directory.

Individual components can be built in-order using the following commands
if changes are required by your build environment

```
make -C forwardmodel/dggev
make -C forwardmodel/spec1d
make -C Reference
make -C InitialPhase/optimizer
make -C Phase/optimizer
'''

## Running

In the tutorial directory, there are a series of numbered bash scripts for running
an example processing of an ambient noise cross-correlation, namely

* 00_create_reference.sh
* 01_create_initial_target_phase.sh
* 02_fit_initial_target_phase.sh
* 03_fit_bessel.sh

Additionally, there is a latex document, tutorial.tex which describes in detail
the inputs and parameters that are used to control the processing of ambient
noise cross correlations with this method.