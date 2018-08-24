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
Eigen::MatrixXcd weights;
double *delays;
double *phases_aligned;

double min_phase_diff_mean = 10*M_PI/180;

//reused buffer
rosjack_data **in_buff;
rosjack_data *out_buff1;
rosjack_data *out_buff2;
Eigen::MatrixXcd in_fft;

double hann(unsigned int buffer_i, unsigned int buffer_size){
    return 0.5 - 0.5*cos(2*PI*buffer_i/(buffer_size-1));
}

void update_weights(bool ini=false){
    int i,j;
    
    double this_dist = 0.0;
    double this_angle = 0.0;
    
    printf("New delays:\n");
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            //assuming first microphone as reference
            delays[i] = 0.0;
            if(ini){
                for(j = 0; j < fft_win; j++){
                    weights(i,j) = 1.0; 
                }
            }
            printf("\t %d -> %f\n",i,delays[i]);
        }else{
            this_dist = array_geometry[i]["dist"];
            this_angle = array_geometry[i]["angle"]-angle;
            if(this_angle>180){
                this_angle -= 360;
            }else if(this_angle<-180){
                this_angle += 360;
            }
            
            delays[i] = this_dist*cos(this_angle*deg2rad)/(-v_sound);
            printf("\t %d -> %f\n",i,delays[i]);
            
            for(j = 0; j < fft_win; j++){
                weights(i,j) = std::exp(-M_I*(double)2*PI*freqs[j]*delays[i]); 
            }
        }
    }
}

double get_overall_phase_diff(int min_i,int *num_i){
    if (min_i < number_of_microphones-1){
        double this_diff = 0;
        for (int i = min_i+1; i < number_of_microphones; i++){
            this_diff += abs(phases_aligned[min_i]-phases_aligned[i]);
            (*num_i)++;
        }
        return this_diff + get_overall_phase_diff(min_i+1,num_i);
    }else{
        return 0;
    }
}

double get_mean(double * data, int data_size){
    double data_sum = 0;
    
    for(int i = 0; i < data_size; i++){
        data_sum += data[i];
    }
    
    return data_sum/(data_size);
}

void apply_weights (rosjack_data **in){
    int i,j;
    double phase_diff_sum;
    int phase_diff_num;
    double phase_diff_mean;
    double mag_mean;
    double pha_mean;
    
    // fft
    for(i = 0; i < number_of_microphones; i++){
        for(j = 0; j < fft_win; j++){
            x_time[j] = in[i][j]*hann_win[j];
        }
        fftw_execute(x_forward);
        for(j = 0; j < fft_win; j++){
            in_fft(i,j) = x_fft[j];
        }
    }
    
    for(j = 0; j < fft_win; j++){
        //applying weights to align phases
        for(i = 0; i < number_of_microphones; i++){
          phases_aligned[i] = arg(conj(weights(i,j))*in_fft(i,j));
        }
        
        //getting the mean phase difference between all microphones
        phase_diff_num = 0;
        phase_diff_sum = get_overall_phase_diff(0,&phase_diff_num);
        phase_diff_mean = phase_diff_sum/(double)phase_diff_num;
        
        if (phase_diff_mean < min_phase_diff_mean){
            //if below threshold, create new frequency data bin
            
            //from mean magnitude
            mag_mean = 0;
            for(i = 0; i < number_of_microphones; i++){
                mag_mean += abs(in_fft(i,j));;
            }
            mag_mean /= number_of_microphones;
            
            //and from mean phase (should be close to 0)
            pha_mean = 0;
            for(i = 0; i < number_of_microphones; i++){
                pha_mean += arg(in_fft(i,j));
            }
            pha_mean /= number_of_microphones;
            
            //TODO: maybe this can be carried out in a quicker way
            y_fft[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
        }else{
            //if not below threshold, nullify data bin
            
            //TODO: this may create artifacts, since it's introducing binary masking
            //      smoothing: no discernable better results
            //      SNR-based masking: ...
            y_fft[j] = 0.0;
        }
    }
    
    // ifft
    fftw_execute(y_inverse);
    
    // preparing output
    for (j = 0; j<fft_win; j++){
        // fftw3 does an unnormalized ifft that requires this normalization
        out_buff2[j] = real(y_time[j])/(double)fft_win;
    }
}

int jack_callback (jack_nframes_t nframes, void *arg){
    //TimeVar t_bef = timeNow();
    int i,j;
    rosjack_data *out_aux;
    
    //Inputing from ROS
    rosjack_data out[nframes];
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        
        for (i = 0; i < number_of_microphones; i++){
            //appending this window to input buffer
            for(j = 0; j < nframes; j++)
                in_buff[i][j+(nframes*3)] = in[i][j];
        }
        
        //applying weights and storing the filter output in out_buff2
        apply_weights(in_buff);
        
        //doing overlap and storing in output
        for(j = 0; j < nframes; j++)
            out[j] = (out_buff1[j+(nframes*2)] + out_buff2[j+nframes])/2;
        
        //storing filtered output to out_buff1 for overlap purposes
        out_aux = out_buff2;
        out_buff2 = out_buff1;
        out_buff1 = out_aux;
        
        //shifting input buffer one window
        for (i = 0; i < number_of_microphones; i++){
            for(j = 0;j < fft_win-nframes; j++)
                in_buff[i][j] = in_buff[i][j+nframes];
        }
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

void theta_roscallback(const std_msgs::Float32::ConstPtr& msg){
    ROS_INFO("Updating weights for angle: %f", msg->data);
    
    angle = msg->data;
    update_weights();
}


int main (int argc, char *argv[]) {
    const char *client_name = "beamform";
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
    out_buff1 = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
    out_buff2 = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
    
    in_buff = (rosjack_data **) malloc (sizeof(rosjack_data*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; i++){
        in_buff[i] = (rosjack_data *) malloc (sizeof(rosjack_data)*fft_win);
    }
    
    freqs = (double *)malloc(sizeof(double)*fft_win);
    for(i = 0; i<fft_win/2;i++){
        freqs[i] = ((double)(i+1)/(double)fft_win)*((double)rosjack_sample_rate);
        freqs[fft_win-1-i] = -((double)(i+1)/(double)fft_win)*((double)rosjack_sample_rate);
    }
    
    delays = (double *) malloc (sizeof(double)*number_of_microphones);
    phases_aligned = (double *) malloc (sizeof(double)*number_of_microphones);
    
    weights.resize(number_of_microphones,fft_win);
    in_fft.resize(number_of_microphones,fft_win);
    
    update_weights(true);
    
    READY = true;
    
    ROS_INFO("Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    ros::spin();
    
    /* apparently ROS requires the use of exit, instead of return for final cleaning up */
    exit(0);
}
