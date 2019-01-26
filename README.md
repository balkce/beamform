# beamform
ROS package that carries out simple 1D beamforming strategies, using JACK as input/output audio server.

It includes a custom ROS topic: JackAudio, defined in the package "jack_msgs". It also uses a custom API for communicating between ROS and JACK referenced as "rosjack", which is not to be confused with https://github.com/balkce/rosjack. Although both serve similar purposes, in this repository, rosjack is built as a library for ROS-JACK inter-operability, not as a ROS package.

Included beamformers:
* das: Delay-and-Sum in frequency domain
* mvdr: Minimum Variance Distortionless Response
* phase: Phase-based binary masking
* gsc: Generalized sidelobe canceller

Two YAML files are required to be configured:
* rosjack_config.yaml: verbosity, type of output (through only JACK, only ROS, or both), if the beamform output should be stored in an audio file, and if the XRUN count should be stored in a an external text file (/home/user/rosjack_xrun_count.txt)
* beamform_config.yaml: used by all beamformers, verbosity, initial direction of interest and microphone positions (required to named as "mic0", "mic1", "mic2", ..., "mic10", etc.).

The direction of interest of all beamformers can be changed on-the-fly by writing to the topic /theta of type std::Float32 (0 is front, -90 is left, 90 is right, 180 is back).

Note for MVDR: it uses only a small portion of the frequencies for speed. It decides which frequencies to use upon a frequency range and basic enery thresholding that can be configured in the mvdr.launch file.

Note for GSC: it uses a dynamic mu that changes depending on the current sample SNR. To facilitate it's configuration, the gsc.launch file includes the values for the starting mu, the maximum mu, the filter size, if a VAD should be used (and its accompanying VAD threshold), and if the behavior of the mu value should be stored in an external text file (/home/user/mu_behavior.txt).

Note for Phase: it process the weight of each frequency bin by considering if the phase is below a threshold; if not, it averages several past magnitudes and multiplies them by a factor to reduce the presence of not-in-phase interferences. To facilitate it's configuration, the phase.launch file includes the values for minimum phase and magnitude factor, and the number of past windows.

Note for audio file: the audio file is a 16-bit WAV file with the sample rate with which the JACK server is configured. If the file path in rosjack_config.yaml is empty, the default path will be used: /home/user/rosjack_write_file.wav.

## Additional utilities included

The following nodes are also included:
* rosjack_read: example of a node that reads from JACK and outputs to the configured output in rosjack_config.yaml (can be the JackAudio ROS topic, JACK, or both).
* rosjack_write: example of a node that writes to JACK from theh JackAudio ROS topic. It uses an intermediate buffer to be robust against inconsistencies in timing, but induces some lag between input and output.

In addition, the rosjack_ref node is also included (with its accompanying launch file) which outputs the reference microphone (mic1) applying the overlap-and-add logistics the beamformers use. These logistics insert a lag which makes it difficult to evaluate the beamformers output online. The rosjack_ref node output should synchronize exactly with the output of the beamformers.

## Dependencies
Packages that can be installed trough apt official repositories:
* libjack-jackd2-dev: JACK development libraries
* libfftw3-dev: a very fast FFT C/C++ implementation
* libeigen3-dev: lightweight linear algebra library
* libsndfile1-dev: library to read and write audio files

