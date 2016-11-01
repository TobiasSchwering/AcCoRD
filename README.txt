
            The AcCoRD Simulator
            (Actor-based Communication via Reaction-Diffusion)
	    with the extension of Cylinders and flow inside them

This document extends the README for AcCoRD v0.5 (public beta, 2016-04-15)

# Introduction

Welcome to the AcCoRD Simulator (Actor-based Communication via Reaction-Diffusion). AcCoRD is a simulation tool for molecular communication. It is a hybrid microscopic/mesoscopic simulator that functions as a generic solver of stochastic reaction-diffusion but which has been developed from the perspective of communications analysis. The focus is to efficiently generate the statistics of molecule observations at some location (i.e., at a receiver).

As of v0.5 (2016-04-15), this document serves as both the official README (https://github.com/adamjgnoel/AcCoRD/blob/master/README.txt) and the Github wiki page (https://github.com/adamjgnoel/AcCoRD/wiki)


## Audience

If you are a communications researcher who is interested in studying molecular communication, then you are the intended audience. However, anyone interested in the multi-scale simulation of reaction-diffusion systems might also find AcCoRD useful.

## Development Status

Development of AcCoRD is on-going and hosted on Github (https://github.com/adamjgnoel/AcCoRD).

The extensions are based on the version 0.5 of the AcCoRD main branch and can be also found on Github (https://github.com/TobiasSchwering/AcCoRD)

# Background

AcCoRD development began in October 2014 with the motivation to develop a reaction-diffusion simulation tool that would be suitable for the communications research community. The goal was to integrate recent advances in stochastic reaction-diffusion models. The overall design of AcCoRD is very much inspired by Smoldyn (http://www.smoldyn.org/), which is also a microscopic reaction-diffusion solver (and which, coincidentally, has also recently added mesoscopic regions). The key differences with Smoldyn relate to AcCoRD's communication focus. AcCoRD's intended use is to capture the dynamics of a communications system, i.e., with a transmitter and receiver. Transmitters are implemented as "active actors", which manually add molecules to the environment. Observations by receivers are recorded as "passive actors". The placement of actors is defined by the user, and the only information that AcCoRD saves is that related to actor behavior. Furthermore, simulations are designed to be run repeatedly (i.e., an arbitrary number of times, but could easily be tens to hundreds of thousands of times), so that accurate receiver statistics can be generated.

AcCoRD is developed in C in order to have precise control over memory management and access to very fast algorithms for random number generation, linked list operations, and related computational tasks. Some utilities have been written (and will be written) in MATLAB to facilitate post-processing.

While a formal publication on the design of AcCoRD is in preparation, the general motivation and some of the implementation ideas were presented in the following conference papers:  
- A. Noel, K.C. Cheung, and R. Schober, On the Statistics of Reaction-Diffusion Simulations for Molecular Communication, in Proc. ACM NANOCOM 2015, Sep. 2015. DOI: http://dx.doi.org/10.1145/2800795.2800821 - the simulations in this paper were completed with proof-of-concept MATLAB code  
- A. Noel, K.C. Cheung, and R. Schober, Multi-Scale Stochastic Simulation for Diffusive Molecular Communication, in Proc. IEEE ICC 2015, pp. 2712--2718, Jun. 2015. DOI: http://dx.doi.org/10.1109/ICC.2015.7248471 - the simulations in this paper were completed with AcCoRD v0.1, which was an early 2D build

The extensions for Cylindric regions and time dependent uniform and laminar flow inside them have been added for the examinations in the master thesis of Tobias Schwering. They are intended to be a fist approximation of blood vessels.

# Feature Summary

The additions extend the main branch AcCoRD by the following features:
- cylindric regions with much integration into the already existing functions. They have to be 'micro' regions
- the possibility to add an axial flow with uniform or laminar profile to cylindrical regions
- time-dependent flow functions, so far linear accel- or decelleration with an initial offset or sinusiodal flow with an optional offset are possible

# Known Issues

see the list of the main branch AcCoRD version 0.5


# Installation

No pre-compiled files are agiven for AcCoRD with these additions. So there is only one intall option:

Build from source code. This is needed if you want to use different compilation
parameters. The src directory contains scripts for Windows, Debian/Ubuntu, and RHEL/CentOS builds (note that the linux build scripts are identical except for the binary filenames in order to distinguish between them). Run a script from
the command line while in the "src" directory and the binary will be placed in the "bin" directory. GCC and standard C libraries are required.
* build_accord_opt_win.bat builds the optimized Windows version with executable accord_win.exe
* build_accord_debug_win.bat builds the debug Windows version with executable accord_win_debug.exe
* build_accord_opt_dub builds the optimized Debian/Ubuntu version with executable accord_dub.out
* build_accord_debug_dub builds the debug Debian/Ubuntu version with executable accord_dub_debug.out
* build_accord_opt_rc builds the optimized RHEL/CentOS version with executable accord_rc.out
* build_accord_debug_rc builds the debug RHEL/CentOS version with executable accord_rc_debug.out

If are compiling on Windows, then minGW is recommended for GCC http://www.mingw.org/

You may need to change the file permissions to execute the build script in Linux (e.g., chmod +x FILENAME).

Example calls for compiling the optimized version:
From Windows command line: build_accord_opt_win.bat
From shell in Debian/Ubuntu: ./build_accord_opt_dub


# Basic Usage

The root source directory includes "fake" config and output files (they include comments so they are not valid JSON files and cannot be used as-is). They are called "HOWTO_DEFINE_CONFIG.txt", "HOWTO_READ_SUMMARY_OUTPUT.txt", and "HOWTO_READ_OUTPUT.txt". The config subdirectory also includes a number of sample valid configuration files for you to start running and experimenting with.

AcCoRD is run from the command line. It takes 2 additional (optional) arguments when called:  
1) Configuration filename. Config file must be defined relative to one of 3 locations. AcCoRD will first search the current directory, then the "config" subdirectory, and finally the "../config/" directory. There are a number of sample configuration files provided to demonstrate AcCoRD functionality.  
2) Seed value for random number generator (will override value specified in config file). You can read more about the meaning of the seed value in HOWTO_DEFINE_CONFIG.txt

If there are no additional arguments, then a default config file is used ("accord_config_sample.txt").
If there is one additional argument, then it must be the configuration filename.

Sample call from Windows command prompt (where both the executable and the configuration file are in the current directory):
accord_win.exe myconfig.txt 2
Sample call from Linux shell while in the AcCoRD "bin" directory while the configuration file is in "../config/":
./accord myconfig.txt 2

AcCoRD will create 2 files. The main output file has the suffix '_SEEDX.txt', where X is the seed number used, and a summary file will have the suffix '_SEEDX_summary.txt'. These 2 files should be kept together, especially if you want to import them into MATLAB.

To import simulation output as a structure in MATLAB, use the accord_import function found in the matlab folder of the source code. Please note that the files in the "JSONLab" subdirectory are required for this function to work. The call is:
[data, config] = accord_import(FILENAME, SEEDRANGE, bWrite)

where
* FILENAME - the output filename to load (including the relative path, if applicable, but excluding the '_SEEDX.txt' suffix)
* SEEDRANGE - a vector specifying the seed values to import (each seed value corresponds to one pair of output files)
* bWrite - switch to write the "data" and "config" structures to a MATLAB mat-file named CONFIG_NAME_out.mat in the "matlab" directory, where CONFIG_NAME is the name of the configuration file that was originally used to run the simulation.
* data - output structure of simulation data
* config - output structure with the parameters of the configuration file used to run the simulation

accord_import will use the summary file to learn the name of the original configuration file and will search for that file in the current directory, then in the "config", "../config/", "../", "../../", and "../../config" directories (and in that order). If a file with a matching configuration filename cannot be found, then "config" will be an empty array.

The "config" structure will contain all of the parameters specified in the original configuration file, using a format similar to that in the configuration file.
The "data" structure will contain the data from both the output and output summary files. The names of the structure members are similar to those used in the output files themselves, so reading the structure should be straightforward.
Examples:
* data.numRepeat -> total number of independent realizations (aggregated from all seed values).
* data.activeID(i) -> index (in the global actor list) of the ith actor in the active actor list. The global list includes both active and passive actors in no particular order.
* data.activeBits{i}(j,k) -> value of the kth bit sent by the ith active actor (indexed from the active actor list) in the jth simulation realization.
* data.passiveRecordMolID{i}(j) -> index (in the global molecule list) of the the jth molecule being observed by the ith recorded passive actor.
* data.passiveRecordCount{i}(j,k,l) -> number of molecules of the kth type observed by the ith recorded passive actor in the lth observation of the jth simulation realization.

A future release of AcCoRD will include more utilities for post-processing simulation results.


# Licensing

Main AcCoRD files are copyright 2016 by Adam Noel under the "New BSD" license. For full details of licensing, including use of third-party code, please see LICENSE.txt

# Credits

Main branch AcCoRD Developed and maintained by Adam Noel (http://www.adamnoel.ca)
Additions Developed and maintained by Tobias Schwering at the Institute for Digital Communications of the University of Erlangen-Nürnberg

Testing: T. Schwering

Supervision and Support of the additions: Prof. R. Schober., Dr. A. Noel, A. Ahmadzadeh, V.Jamali