/**
 * Agent that outputs the first microphone of the array in the same delayed manner as the rest of the frequency-domain beamformers do when using the overlap-and-add logistics
 */

#include "rosjack.h"
#include "util.h"

#include <map>
#include <vector>
#include <cmath>

// Include FFTW header
#include <complex>
#include <fftw3.h>

bool READY = false;
std::complex<double> * x;

void apply_weights (jack_ringbuffer_t *in, rosjack_data *out, int mic){
    int j;
    
    overlap_and_add_prepare_input(in, x);
    
    for (j = 0; j<fft_win; j++){
        out[j] = real(x[j]);
        
        //applying wola to avoid discontinuities in the time domain
        out[j] *= hann_win[j];
    }
}

int jack_callback (jack_nframes_t nframes, void *arg){
    TimeVar t_bef = timeNow();
    int i,j;
    
    //Inputing from ROS
    rosjack_data **out;  //needs to be dynamically allocated to be used by do_overlap_bymic
    out = (rosjack_data **) malloc (sizeof(rosjack_data*)*1);
    out[0] = (rosjack_data *) calloc (nframes,sizeof(rosjack_data));
    
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        do_overlap_bymic(in, out, nframes, apply_weights);
    }else{
        for (i = 0; i < nframes; i++){
            out[0][i] = 0.0;
        }
    }
    
    //Outputing to ROS
    output_to_rosjack (out[0], nframes, output_type);
    
    //std::cout << "Callback took: " << duration(timeNow()-t_bef)/1000000.0 << " ms.\n";
    return 0;
}

int main (int argc, char *argv[]) {
    const char *client_name_local = "rosjack_ref";
    int i;
    /* ROS initialization*/
    ros::init(argc, argv, client_name_local);
    ros::NodeHandle n;
    handle_params(&n);
    
    number_of_microphones = 1;
    
    /* create JACK agent */
    if(rosjack_create (ROSJACK_READ, &n, "jackaudio_ref", client_name_local, number_of_microphones, jack_callback)){
        ROS_ERROR("JACK agent could not be created.\n");
        ros::shutdown();
        exit(1);
    }
    
    std::cout << "Pre-allocating space for internal buffers." << std::endl;
    prepare_overlap_and_add_bymic(); //fft_win is assinged here
    
    x = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    
    READY = true;
    
    ROS_INFO("Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    ros::spin();
    
    /* apparently ROS requires the use of exit, instead of return for final cleaning up */
    exit(0);
}
