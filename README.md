# beamform
ROS package that carries out simple 1D beamforming strategies, using JACK as input/output audio server.

It includes a custom ROS topic: JackAudio, defined in the package "jack_msgs". It also uses a custom API for communicating between ROS and JACK referenced as "rosjack", which is not to be confused with https://github.com/balkce/rosjack. Although both serve similar purposes, in this repository, rosjack is built as a library for ROS-JACK inter-operability, not as a ROS package.

Included beamformers:
* das: Delay-and-Sum in frequency domain
* mvdr: Minimum Variance Distortionless Response
* gsc: Generalized sidelobe canceller
* lcmv: Linearly Constrained Minimum Variance
* gss: Geometric Source Separation
* phase: Phase-based binary masking
* phasempf: Phase-based binary masking with post-filter

Two YAML files are required to be configured:
* rosjack_config.yaml: verbosity, type of output (through only JACK, only ROS, or both), if the beamform output should be stored in an audio file, if the XRUN count should be stored in a an external text file (/home/user/rosjack_xrun_count.txt), and the sample rate for the ROS output (if audio file is being stored, the ROS sample rate is used, if not, the JACK sample rate is used).
* beamform_config.yaml: used by all beamformers, verbosity, initial direction of interest and microphone positions (required to named as "mic0", "mic1", "mic2", ..., "mic10", etc.).

The direction of interest of all beamformers can be changed on-the-fly by writing to the topic /theta of type std::Float32 (0 is front, -90 is left, 90 is right, 180 is back).

Note for MVDR: it uses only a small portion of the frequencies for speed. It decides which frequencies to use upon a frequency range and basic enery thresholding that can be configured in the mvdr.launch file.

Note for GSC: it uses a dynamic mu that changes depending on the current sample SNR. To facilitate it's configuration, the gsc.launch file includes the values for the starting mu, the maximum mu, the filter size, if a VAD should be used (and its accompanying VAD threshold), and if the behavior of the mu value should be stored in an external text file (/home/user/mu_behavior.txt).

Note for LCMV: similar to MVDR, it uses only a small portion of the frequencies for speed. It also considers the positions of the interferences for nullifying effect. For this, beamform_config.yaml also stores the initial direction of interferences, and LCMV ignores the interference that have an absolute value greater than 180 and the ones that follow it in the interference list. Similar to the /theta topic, LCMV also listens for the /theta_interference topic, that uses a custom InterfTheta message with the following structure {id,angle}. If the interference id is outside the range [1,interference number] it will add another interference to the list. If the new angle of a current interference is too close to a current interference, it is eliminated.

Note for GSS: similar to MVDR and LCMV, it uses only a small portion of the frequencies for speed (although it is fast enough such that the magnitude threshold can be considerably lower than MVDR and LCMV). Also similar to LCMV, it also considers the positions of the interferences, however it uses them, as well as the position of the source of interest, for its inner processes. The workings of how interferences are inserted and removed are the same as in LCMV. Although in its original form GSS is able to separate all the sources (of interest and interferences), in this implementation it only provides the one of interest since it would require the reworking of the whole ROS package to consider more than one output.

Note for Phase: it process the weight of each frequency bin by considering if the inter-microphone phase difference is below a threshold; if not, it multiplies the reference microphone magnitude by a factor to reduce the presence of not-in-phase interferences. To facilitate it's configuration, the phase.launch file includes the values for minimum phase and magnitude factor, and the size of its output smoothing filter.

Note for PhaseMPF: carries out the Phase beamformer and an anti-Phase beamfomer (by using the negative mask). Both outputs are then inputed into a bi-channel post-filter, which is a variation based on the multi-channel post-filter of Valin et al (2007), the source code of which can be found in https://github.com/introlab/manyears/blob/master/manyears-C/dsplib/Separation/postfilter.c. In turn, this post-filter uses the Minima Controlled Recursive Averaging (MCRA) noise estimator of Cohen and Berdugo (2001). To facilitate its configuration the phasempf.launch includes all the parameters of the Phase beamformer and the modified post-filter, including its MCRA module.

Note for audio file: the audio file is a 16-bit WAV file with the sample rate with which the JACK server is configured. If the file path in rosjack_config.yaml is empty, the default path will be used: /home/user/rosjack_write_file.wav.

## Additional utilities included

The following nodes are also included:
* rosjack_read: example of a node that reads from JACK and outputs to the configured output in rosjack_config.yaml (can be the JackAudio ROS topic, JACK, or both).
* rosjack_write: example of a node that writes to JACK from the JackAudio ROS topic. It uses an intermediate buffer to be robust against inconsistencies in timing, but induces some lag between input and output.
* rosjack_ref: (with accompanying launch file) outputs the reference microphone (mic1) applying the overlap-and-add logistics the beamformers use. These logistics insert a lag which makes it difficult to evaluate the beamformers output online. The rosjack_ref node output should synchronize exactly with the output of the beamformers.
* mcra: (with accompanying launch file) carries out the Minima Controlled Recursive Averaging (MCRA) noise estimator with the mic1 input, it removes the noise estimation from mic1 and outputs the result.

## Dependencies
Packages that can be installed trough apt official repositories:
* libjack-jackd2-dev: JACK development libraries
* libfftw3-dev: a very fast FFT C/C++ implementation
* libeigen3-dev: lightweight linear algebra library
* libsndfile1-dev: library to read and write audio files
* libsamplerate0-dev: library for sample rate conversion

## References
Some useful references:
* J. Valin, S. Yamamoto, J. Rouat, F. Michaud, K. Nakadai and H. G. Okuno, "Robust Recognition of Simultaneous Speech by a Mobile Robot," in IEEE Transactions on Robotics, vol. 23, no. 4, pp. 742-752, 2007.
* I. Cohen, and B. Berdugo,"Speech enhancement for non-stationary noise environments," in Signal Processing, vol. 81, no. 11, pp. 2403-2418, 2001.

