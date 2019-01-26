/**
 * ROS agent that reads from microphone and outputs to ROS topic
 */

#include "rosjack.h"
#include "util.h"

#include <map>
#include <vector>
#include <cmath>

// Include FFTW header
#include <complex>
#include <fftw3.h>

// Eigen include
#include <Eigen/Eigen>

bool READY = false;

std::complex<double> *x_fft, *x_time, *y_fft, *y_time;
fftw_plan x_forward, y_inverse;

double *hann_win;
double *freqs;
unsigned int fft_win;

//reused buffers
rosjack_data *in_buff;
rosjack_data *out_buff1;
rosjack_data *out_buff2;

double hann(unsigned int buffer_i, unsigned int buffer_size){
    return 0.5 - 0.5*cos(2*PI*buffer_i/(buffer_size-1));
}

int jack_callback (jack_nframes_t nframes, void *arg){
    TimeVar t_bef = timeNow();
    int i,j;
    
    //Inputing from ROS
    rosjack_data out[nframes];
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        
        //appending this window to input buffer
        for(j = 0; j < nframes; j++)
            in_buff[j+(nframes*3)] = in[0][j];
        
        //applying hann window
        memcpy(out_buff2, in_buff,sizeof(rosjack_data)*fft_win);
        for(j = 0; j < fft_win; j++){
            out_buff2[j] *= hann_win[j];
        }

        //doing overlap and storing in output
        for(j = 0; j < nframes; j++)
            out[j] = out_buff1[j+(nframes*2)] + out_buff2[j+nframes];
        
        //storing filtered output to out_buff1 for overlap purposes
        memcpy(out_buff1, out_buff2,sizeof(rosjack_data)*fft_win);
        
        //shifting input buffer one window
        for(j = 0;j < fft_win-nframes; j++)
            in_buff[j] = in_buff[j+nframes];
    }else{
        for (i = 0; i < nframes; i++){
            out[i] = 0.0;
        }
    }
    
    //Outputing to ROS
    output_to_rosjack (out, nframes, output_type);
    
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
    fft_win = rosjack_window_size*4;
    
    x_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    x_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    
    x_forward = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(x_time), reinterpret_cast<fftw_complex*>(x_fft), FFTW_FORWARD, FFTW_MEASURE);
    y_inverse = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(y_fft), reinterpret_cast<fftw_complex*>(y_time), FFTW_BACKWARD, FFTW_MEASURE);
    
    hann_win = (double *) malloc(sizeof(double) * fft_win);
    for (i = 0; i < fft_win; i++){
        hann_win[i] = hann(i, fft_win);
    }
    out_buff1 = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
    out_buff2 = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
    
    in_buff = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
    
    READY = true;
    
    ROS_INFO("Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    ros::spin();
    
    /* apparently ROS requires the use of exit, instead of return for final cleaning up */
    exit(0);
}
