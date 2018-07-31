# beamform
ROS package that carries out simple 1D beamforming strategies, using JACK as input/output audio server.

It includes a custom ROS topic: JackAudio, defined in the package "jack_msgs". It also uses a custom API for communicating between ROS and JACK referenced as "rosjack", which is not to be confused with https://github.com/balkce/rosjack. Although both serve similar purposes, in this repository, rosjack is built as a library for ROS-JACK inter-operability, not as a ROS package.

Included beamformers:
* das_time: Delay-and-Sum in time domain
* das: Delay-and-Sum in frequency domain
* mvdr: Minimum Variance Distortionless Response
* phase: Phase-based binary masking

Two YAML files are required to be configured:
* rosjack_config.yaml: verbosity and type of output (through only JACK, only ROS, or both)
* beamform_config.yaml: used by all beamformers, verbosity, initial direction of interest and microphone positions (required to named as "mic0", "mic1", "mic2", ..., "mic10", etc.).

The direction of interest of all beamformers can be changed on-the-fly by writing to the topic /theta of type std::Float32 (0 is front, -90 is left, 90 is right, 180 is back).

Note for MVDR: it uses only a small portion of the frequencies for speed. It decides which frequencies to use upon a hardcoded frequency range and basic enery thresholding.

## Dependencies
Packages that can be installed trough apt official repositories:
* libjack-jackd2-dev: JACK development libraries
* libfftw3-dev: a very fast FFT C/C++ implementation
* libeigen3-dev: lightweight linear algebra library

