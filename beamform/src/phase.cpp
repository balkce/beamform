/**
 * Beamform that carries out phase-based frequency masking by simple thresholding.
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
Eigen::MatrixXcd weights;

double min_mag = 0.0001;
double min_phase = 10;
double min_phase_diff_mean = 10*M_PI/180;
int smooth_size = 10;
double past_out = 0.0;

//reused buffers
double *phases_aligned;
double *past_samples;
Eigen::MatrixXcd in_fft;

void shift_data(double data, double *buf, int size){
  for(int i = 1; i < size; i++){
    buf[i-1] = buf[i];
  }
  buf[size-1] = data;
}

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
}

double get_overall_phase_diff(int min_i,int *num_i){
    if (min_i < number_of_microphones-1){
        double this_diff = 0;
        double this_diff_raw;
        for (int i = min_i+1; i < number_of_microphones; i++){
            this_diff_raw = abs(phases_aligned[min_i]-phases_aligned[i]);
            if (this_diff_raw > M_PI)
                this_diff_raw = 2*M_PI - this_diff_raw;
            this_diff += this_diff_raw;
            (*num_i)++;
        }
        return this_diff + get_overall_phase_diff(min_i+1,num_i);
    }else{
        return 0;
    }
}

double get_mean(double * data, int data_size){
    double data_sum = 0.0;
    
    for(int i = 0; i < data_size; i++){
        data_sum += data[i];
    }
    
    return data_sum/(double)data_size;
}

void apply_weights (rosjack_data **in, rosjack_data *out){
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
    
    y_fft[0] = in_fft(0,0);
    for(j = 1; j < fft_win; j++){
        //applying weights to align phases
        for(i = 0; i < number_of_microphones; i++){
          phases_aligned[i] = arg(conj(weights(i,j))*in_fft(i,j));
        }
        
        //getting the mean phase difference between all microphones
        phase_diff_num = 0;
        phase_diff_sum = get_overall_phase_diff(0,&phase_diff_num);
        phase_diff_mean = phase_diff_sum/(double)phase_diff_num;
        
        //creating new frequency data bin from mean magnitude
        mag_mean = 0;
        for(i = 0; i < number_of_microphones; i++){
            mag_mean += abs(in_fft(i,j));
        }
        mag_mean /= number_of_microphones;
        
        //and from the phase of the reference microphone
        pha_mean = arg(in_fft(0,j));
        
        
        if (phase_diff_mean < min_phase_diff_mean){
            y_fft[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
        }else{
            mag_mean *= min_mag;
            y_fft[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
        }
        
        
        //mag_mean *= 1 / (1+(15.0)*phase_diff_mean); // 15.0 here is a constant that is linked to the minimum phase difference
        //y_fft[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
    }
    
    // ifft
    fftw_execute(y_inverse);
    
    // preparing output
    for (j = 0; j<fft_win; j++){
        // fftw3 does an unnormalized ifft that requires this normalization
        out[j] = real(y_time[j])/(double)fft_win;
    }
}

int jack_callback (jack_nframes_t nframes, void *arg){
    //TimeVar t_bef = timeNow();
    int i,j;
    
    //Inputing from ROS
    rosjack_data out[nframes];
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        do_overlap(in, out, nframes, apply_weights);
    
        // this type of masking requires to smooth the output since it
        // inserts many discontinuities in the time domain
        for (j = 0; j<nframes; j++){
            shift_data(out[j], past_samples,smooth_size);
            out[j] = get_mean(past_samples,smooth_size);
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

void phase_handle_params(ros::NodeHandle *n){
    std::string node_name = ros::this_node::getName();
    std::cout << "Phase ROS parameters: " << std::endl;
    
    if ((*n).getParam(node_name+"/min_phase",min_phase)){
        ROS_INFO("Min Phase Threshold: %f",min_phase);
    }else{
        min_phase = 10.0;
        ROS_WARN("Min Phase Threshold argument not found in ROS param server, using default value (%f).",min_phase);
    }
    min_phase_diff_mean = min_phase*M_PI/180;
    
    if ((*n).getParam(node_name+"/min_mag",min_mag)){
        ROS_INFO("Min Mag Threshold: %f",min_mag);
    }else{
        min_mag = 10.0;
        ROS_WARN("Min Mag Threshold argument not found in ROS param server, using default value (%f).",min_mag);
    }
    
    if ((*n).getParam(node_name+"/smooth_size",smooth_size)){
        ROS_INFO("Smooth Filter Size: %d",smooth_size);
        if(smooth_size < 1){
          smooth_size = 20;
          ROS_WARN("Invalid Smooth Filter Size argument, using default value (%d).",smooth_size);
        }
    }else{
        smooth_size = 20;
        ROS_WARN("Smooth Filter Size argument not found in ROS param server, using default value (%d).",smooth_size);
    }
}

int main (int argc, char *argv[]) {
    int i;
    /* ROS initialization*/
    ros::init(argc, argv, client_name);
    ros::NodeHandle n;
    handle_params(&n);
    phase_handle_params(&n);
    
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
    
    past_samples = (double *) malloc (sizeof(double)*smooth_size);
    for (i = 0; i < smooth_size; i++){
        past_samples[i] = 0.0;
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
