# jarvis

Software to carry out Principal Component Analysis (PCA) of experimental (or simulated) data.

This software was originally developed by Jorge Quintanilla as part of two collaborations involving Stuart Gibson, Robert Twyman, Dylan Barker, Gunnar Moller and Tymoteusz Tula. For licensing information, see COPYING. You can cite this code using the following DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4266743.svg)](https://doi.org/10.5281/zenodo.4266743)

This software is optimised for Octave but can be easily adapted for use with Matlab.

REQUIREMENTS

A working installation of Octave (the software has been tested on Octave 3.4.3; the source code has comments explaining how it can be adapted for use on Matlab).

INSTALLATION

No installation required - just place the file functions_PCA_current.m in your working directory.

USAGE

Start Octave.

Type the command

 source "functions_PCA_current.m"

This will make available the following Octave functions:

PCA
gen_A
gen_image_matrix
recons
scores

The usage of each of these functions is documented thorugh comments within its source code.
