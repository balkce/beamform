/**
 * ROS agent that reads from microphone and outputs to ROS topic
 */

#include "rosjack.h"
#include "util.h"

// Include FFTW header
#include <complex>
#include <fftw3.h>

bool READY = false;

std::complex<double> *x_fft, *x_time, *y_fft, *y_time;
fftw_plan x_forward, y_inverse;

double *hann_win;
double *freqs;
unsigned int fft_win;
std::complex<double> **weights;
double *delays;

rosjack_data **in_overlap;
rosjack_data **in_overlap2;

double hann(unsigned int buffer_i, unsigned int buffer_size){
    return 0.5 - 0.5*cos(2*PI*buffer_i/(buffer_size-1));
}

void update_weights(bool ini=false){
    calculate_delays(delays);
    
    int i,j;
    
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            if(ini){
                for(j = 0; j < fft_win; j++){
                    weights[i][j] = 1.0; 
                }
            }
        }else{
            for(j = 0; j < fft_win; j++){
                weights[i][j] = std::exp(-M_I*(double)2*PI*freqs[j]*delays[i]); 
            }
        }
    }
}

rosjack_data * apply_weights (rosjack_data *in, int mic){
    int i;

    // fft
    for(i = 0; i < fft_win; i++){
        x_time[i] = in[i]*hann_win[i];
    }
    fftw_execute(x_forward);

    // applying weights per frequency
    for(i = 0; i < fft_win; i++){
        x_fft[i] *= conj(weights[mic][i]); 
    }

    // ifft
    for(i = 0; i < fft_win; i++){
        y_fft[i] = x_fft[i];
    }
    fftw_execute(y_inverse);

    // preparing output
    rosjack_data * out = (rosjack_data *)malloc(sizeof(rosjack_data)*fft_win);
    for (i = 0; i<fft_win; i++){
        // fftw3 does an unnormalized ifft that requires this normalization
        out[i] = real(y_time[i])/(double)fft_win;
    }

    return out;
}

int jack_callback (jack_nframes_t nframes, void *arg){
    int i,j;
    
    //Inputing from ROS
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        rosjack_data aux[number_of_microphones][fft_win];
        rosjack_data out[nframes];
        
        
        for(j = 0;j < nframes; j++)
            out[j] = 0.0;
            
        for (i = 0; i < number_of_microphones; i++){
            //adding this window to overlap buffer2
            for(j = 0; j < nframes; j++)
                in_overlap2[i][j+(nframes*3)] = in[i][j];
            
            //applying weights
            memcpy(aux[i], apply_weights(in_overlap2[i],i),sizeof(rosjack_data)*fft_win);
            
            //doing overlap
            for(j = 0; j < nframes; j++)
                out[j] += (in_overlap[i][j+(nframes*2)] + aux[i][j+nframes])/number_of_microphones;
        }
        
        //Outputing to ROS
        output_to_rosjack (out, nframes, output_type);
        
        //copying aux to overlap buffer
        for (i = 0; i < number_of_microphones; i++){
            memcpy(in_overlap[i], aux[i],sizeof(rosjack_data)*fft_win);
        }
        //shifting in_overlap buffer
        for (i = 0; i < number_of_microphones; i++){
            for(j = 0;j < fft_win-nframes; j++)
                in_overlap2[i][j] = in_overlap2[i][j+nframes];
        }
    }else{
        rosjack_data zeros[1024];
        for (i = 0; i < nframes; i++){
            zeros[i] = 0.0;
        }
        output_to_rosjack (zeros, nframes, output_type);
    }
    return 0;
}

void theta_roscallback(const std_msgs::Float32::ConstPtr& msg){
    ROS_INFO("Updating weights for angle: %f", msg->data);
    
    angle = msg->data;
    update_weights();
}


int main (int argc, char *argv[]) {
    int i;
    /* ROS initialization*/
    ros::init(argc, argv, client_name);
    ros::NodeHandle n;
    handle_params(&n);
    
    ros::Subscriber theta_subscriber = n.subscribe("theta", 1000, theta_roscallback);
    
    
    /* create JACK agent */
    if(rosjack_create (ROSJACK_READ, &n, "jackaudio", client_name, number_of_microphones, jack_callback)){
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
    
    in_overlap = (rosjack_data **) malloc (sizeof(rosjack_data*)*number_of_microphones);
    in_overlap2 = (rosjack_data **) malloc (sizeof(rosjack_data*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; i++){
        in_overlap[i] = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
        in_overlap2[i] = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
    }
    
    freqs = (double *)malloc(sizeof(double)*fft_win);
    calculate_frequency_vector(freqs,fft_win);
    
    weights = (std::complex<double> **) malloc (sizeof(std::complex<double>*)*number_of_microphones);
    for(i = 0; i<number_of_microphones;i++){
        weights[i] = (std::complex<double> *) malloc(sizeof(std::complex<double>) * fft_win);
    }
    
    delays = (double *) malloc (sizeof(double)*number_of_microphones);
    
    update_weights(true);
    
    READY = true;
    
    ROS_INFO("Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    ros::spin();
    
    /* apparently ROS requires the use of exit, instead of return for final cleaning up */
    exit(0);
}
