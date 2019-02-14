/**
 * Mininimum-variance distortionless response, with frequency magnitude
 * thresholding so that it can be carried out online.
 */

#include "rosjack.h"
#include "util.h"

// Include FFTW header
#include <complex>
#include <fftw3.h>

// Eigen include
#include <Eigen/Eigen>

bool READY = false;

std::complex<double> *x_fft, *x_time, *y_fft, *y_time;
fftw_plan x_forward, y_inverse;

double *freqs;
double *delays;
Eigen::MatrixXcd weights(1,1);
Eigen::MatrixXcd weights_h(1,1);

//MVDR config
unsigned int past_windows = 10;
double freq_mag_threshold = 1.5;
double freq_max = 4000;
double freq_min = 400;
double out_amp = 4.5;

//reused buffers
Eigen::MatrixXcd in_fft(1,1);
Eigen::MatrixXcd R(1,1);
Eigen::MatrixXcd invR(1,1);
Eigen::MatrixXcd whiteR(1,1);
Eigen::MatrixXcd MVDR_w(1,1);
std::vector<Eigen::MatrixXcd> past_ffts;

void update_weights(bool ini=false){
    calculate_delays(delays);
    
    int i,j;
    
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            if(ini){
                for(j = 0; j < fft_win; j++){
                    weights(i,j) = 1.0; 
                }
            }
        }else{
            for(j = 0; j < fft_win; j++){
                weights(i,j) = std::exp(-M_I*(double)2*PI*freqs[j]*delays[i]); 
            }
        }
    }
    weights_h = weights.adjoint();
}

void apply_weights (rosjack_data **in, rosjack_data *out){
    int i,j;
    double this_freq, this_mag;
    
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
    
    //applying weights
    for(j = 0; j < fft_win; j++){
        this_freq = abs(freqs[j]);
        this_mag = 0.0;
        for(i = 0; i < number_of_microphones; i++)
            this_mag += abs(in_fft(i,j));
        this_mag /= number_of_microphones*fft_win; //normalize the average against window length
        
        if(this_freq >= freq_min && this_freq <= freq_max){
            if(this_mag > freq_mag_threshold){
                //calculating sample covariance for this fft bin
                R = (past_ffts[j]*past_ffts[j].adjoint()).cwiseProduct(whiteR);
                invR = R.inverse();
                
                //calculating MVDR optimal weights for this fft bin
                MVDR_w = (invR*weights.col(j))/((weights_h.row(j)*invR*weights.col(j)));
                
                //applying weigths
                y_fft[j] = (MVDR_w.adjoint()*in_fft.col(j))(0,0);
            }else{
                y_fft[j] = in_fft(0,j)*0.01;
            }
            
            //shifting past_ffts and appending this window
            past_ffts[j].block(0,0,number_of_microphones,past_windows-1) = past_ffts[j].block(0,1,number_of_microphones,past_windows-1);
            past_ffts[j].col(past_windows-1) = in_fft.col(j);
        }else{
            y_fft[j] = 0.0;
        }
    }
    
    // ifft
    fftw_execute(y_inverse);
    
    // preparing output
    for (j = 0; j<fft_win; j++){
        // fftw3 does an unnormalized ifft that requires this normalization
        out[j] = (real(y_time[j])/(double)fft_win)*out_amp;
    }
}

int jack_callback (jack_nframes_t nframes, void *arg){
    TimeVar t_bef = timeNow();
    int i,j;
    
    //Inputing from ROS
    rosjack_data out[nframes];
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        do_overlap(in, out, nframes, apply_weights);
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

void mvdr_handle_params(ros::NodeHandle *n){
    std::string node_name = ros::this_node::getName();
    std::cout << "MVDR ROS parameters: " << std::endl;
    
    double past_windows_tmp;
    if ((*n).getParam(node_name+"/past_windows",past_windows_tmp)){
        past_windows = (int)past_windows_tmp;
        ROS_INFO("Number of Windows for Covariance Calculation: %d",past_windows);
    }else{
        past_windows = 10;
        ROS_WARN("Number of Windows for Covariance Calculation argument not found in ROS param server, using default value (%d).",past_windows);
    }
    
    if ((*n).getParam(node_name+"/freq_mag_threshold",freq_mag_threshold)){
        ROS_INFO("Frequency Magnitude Threshold: %f",freq_mag_threshold);
    }else{
        freq_mag_threshold = 1.5;
        ROS_WARN("Frequency Magnitude Threshold argument not found in ROS param server, using default value (%f).",freq_mag_threshold);
    }
    
    if ((*n).getParam(node_name+"/freq_max",freq_max)){
        ROS_INFO("Max Frequency: %f",freq_max);
    }else{
        freq_max = 4000;
        ROS_WARN("Max Frequency argument not found in ROS param server, using default value (%f).",freq_max);
    }
    
    if ((*n).getParam(node_name+"/freq_min",freq_min)){
        ROS_INFO("Min Frequency: %f",freq_min);
    }else{
        freq_min = 400;
        ROS_WARN("Min Frequency argument not found in ROS param server, using default value (%f).",freq_min);
    }
    
    if ((*n).getParam(node_name+"/out_amp",out_amp)){
        ROS_INFO("Output Amplification: %f",out_amp);
    }else{
        out_amp = 4.5;
        ROS_WARN("Output Amplification argument not found in ROS param server, using default value (%f).",out_amp);
    }
    
}


int main (int argc, char *argv[]) {
    int i;
    /* ROS initialization*/
    ros::init(argc, argv, client_name);
    ros::NodeHandle n;
    handle_params(&n);
    mvdr_handle_params(&n);
    
    ros::Subscriber theta_subscriber = n.subscribe("theta", 1000, theta_roscallback);
    
    
    /* create JACK agent */
    if(rosjack_create (ROSJACK_READ, &n, "jackaudio", client_name, number_of_microphones, jack_callback)){
        ROS_ERROR("JACK agent could not be created.\n");
        ros::shutdown();
        exit(1);
    }
    
    std::cout << "Pre-allocating space for internal buffers." << std::endl;
    prepare_overlap_and_add(); //fft_win is assinged here
    
    x_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    x_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    
    x_forward = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(x_time), reinterpret_cast<fftw_complex*>(x_fft), FFTW_FORWARD, FFTW_MEASURE);
    y_inverse = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(y_fft), reinterpret_cast<fftw_complex*>(y_time), FFTW_BACKWARD, FFTW_MEASURE);
    
    freqs = (double *)malloc(sizeof(double)*fft_win);
    calculate_frequency_vector(freqs,fft_win);
    
    delays = (double *) malloc (sizeof(double)*number_of_microphones);
    
    weights.resize(number_of_microphones,fft_win);
    weights_h.resize(fft_win,number_of_microphones);
    in_fft.resize(number_of_microphones,fft_win);
    
    for(i = 0; i<fft_win;i++){
        Eigen::MatrixXcd this_past_fs(number_of_microphones,past_windows);
        this_past_fs.setZero();
        past_ffts.push_back(this_past_fs);
    }
    
    MVDR_w.resize(number_of_microphones,1);
    
    R.resize(number_of_microphones,number_of_microphones);
    invR.resize(number_of_microphones,number_of_microphones);
    
    whiteR.resize(number_of_microphones,number_of_microphones);
    whiteR.setOnes();
    for(i = 0; i<number_of_microphones;i++){
        whiteR(i,i) = 1.001;
    }
    
    update_weights(true);
    
    READY = true;
    
    ROS_INFO("Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    ros::spin();
    
    /* apparently ROS requires the use of exit, instead of return for final cleaning up */
    exit(0);
}
