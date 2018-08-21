# beamform
ROS package that carries out simple 1D beamforming strategies, using JACK as input/output audio server.

It includes a custom ROS topic: JackAudio, defined in the package "jack_msgs". It also uses a custom API for communicating between ROS and JACK referenced as "rosjack", which is not to be confused with https://github.com/balkce/rosjack. Although both serve similar purposes, in this repository, rosjack is built as a library for ROS-JACK inter-operability, not as a ROS package.

Included beamformers:
* das_time: Delay-and-Sum in time domain
* das: Delay-and-Sum in frequency domain
* mvdr: Minimum Variance Distortionless Response
* phase: Phase-based binary masking
* gsc: Generalized sidelobe canceller

Two YAML files are required to be configured:
* rosjack_config.yaml: verbosity, type of output (through only JACK, only ROS, or both), if the beamform output should be stored in an audio file, and if the XRUN count should be stored in a an external text file (/home/user/rosjack_xrun_count.txt)
* beamform_config.yaml: used by all beamformers, verbosity, initial direction of interest and microphone positions (required to named as "mic0", "mic1", "mic2", ..., "mic10", etc.).

The direction of interest of all beamformers can be changed on-the-fly by writing to the topic /theta of type std::Float32 (0 is front, -90 is left, 90 is right, 180 is back).

Note for MVDR: it uses only a small portion of the frequencies for speed. It decides which frequencies to use upon a hardcoded frequency range and basic enery thresholding.

Note for GSC: it uses a dynamic mu that changes depending on the current sample SNR. To facilitate it's configuration, the YAML file gsc_config.yaml is provided to establish the starting mu, the maximum mu, the filter size, if a VAD should be used (and its accompanying VAD threshold), and if the behavior of the mu value should be stored in an external text file (/home/user/mu_behavior.txt).

Note for audio file: the audio file is a 16-bit WAV file with the sample rate with which the JACK server is configured. If the file path in rosjack_config.yaml is empty, the default path will be used: /home/user/rosjack_write_file.wav.

## Dependencies
Packages that can be installed trough apt official repositories:
* libjack-jackd2-dev: JACK development libraries
* libfftw3-dev: a very fast FFT C/C++ implementation
* libeigen3-dev: lightweight linear algebra library
* libsndfile1-dev: library to read and write audio files

