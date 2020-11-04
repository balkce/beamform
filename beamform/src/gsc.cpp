/**
 * Generalized sidelobe canceller, with a dynamic adaptation rate (mu).
 */

#include "rosjack.h"
#include "util.h"

// Include FFTW header
#include <complex>
#include <fftw3.h>

bool READY = false;

std::complex<double> *x_fft, *x_time, *y_fft, *y_time;
fftw_plan x_forward, y_inverse;

bool use_vad;
double vad_threshold;
double mu0;
double mu_max;
int filter_size;
bool write_mu;
char mu_file_path[FILENAME_MAX];
FILE *mu_file;
rosjack_data **block_matrix;
rosjack_data **filter;
rosjack_data *last_outputs;
rosjack_data last_avg_mu = 0.0;

double *freqs;
double *delays;
std::complex<double> **weights;

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

void apply_weights (jack_ringbuffer_t *in, rosjack_data *out, int mic){
    int i;
    
    // fft
    overlap_and_add_prepare_input(in, x_time);
    fftw_execute(x_forward);
    
    // applying weights per frequency
    for(i = 0; i < fft_win; i++){
        x_fft[i] *= conj(weights[mic][i]);
        //x_fft[i] /= number_of_microphones;
    }
    
    // ifft
    for(i = 0; i < fft_win; i++){
        y_fft[i] = x_fft[i];
    }
    fftw_execute(y_inverse);
    
    // preparing output
    overlap_and_add_prepare_output(y_time,out);
}

void shift_data(rosjack_data data, rosjack_data *buf, int size){
  for(int i = 1; i < size; i++){
    buf[i-1] = buf[i];
  }
  buf[size-1] = data;
}

rosjack_data calculate_power(rosjack_data *buf, int size){
  rosjack_data p = 0.0;
  for(int i = 0; i < size; i++){
    p += buf[i]*buf[i];
  }
  p /= (rosjack_data)size;
  return sqrt(p);
}

int jack_callback (jack_nframes_t nframes, void *arg){
    int i,j,k;
    
    //Inputing from ROS
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        rosjack_data aux[number_of_microphones][fft_win];
        rosjack_data das_out;
        rosjack_data block_out;
        rosjack_data last_out_power;
        rosjack_data block_power;
        rosjack_data out[nframes];
        rosjack_data this_mu;
        rosjack_data avg_mu = 0;
        
        rosjack_data **overlap_out; //needs to be dynamically allocated to be used by do_overlap_bymic
        overlap_out = (rosjack_data **) malloc (sizeof(rosjack_data*)*number_of_microphones);
        for (i = 0; i < number_of_microphones; ++i){
            overlap_out[i] = (rosjack_data *) calloc (nframes,sizeof(rosjack_data));
        }
        
        for(j = 0;j < nframes; ++j)
          out[j] = 0.0;
        
        do_overlap_bymic(in, overlap_out, nframes, apply_weights);
        
        //doing GSC for each sample
        for(j = 0; j < nframes; ++j){
          //doing the upper beamform
          das_out = 0.0;
          for (i = 0; i < number_of_microphones; ++i)
            das_out += overlap_out[i][j];
          das_out /= number_of_microphones;
          
          out[j] = das_out;
          for (i = 0; i < number_of_microphones-1; ++i){
            //blocking matrix
            shift_data(overlap_out[i+1][j]-overlap_out[i][j],block_matrix[i],filter_size);
            
            //applying filters
            block_out = 0.0;
            for(k = 0; k < filter_size; ++k){
              block_out += filter[i][k]*block_matrix[i][k];
            }
            
            //applying to output
            out[j] -= block_out;
          }
          
          //updating last outputs and their power
          shift_data(out[j],last_outputs,filter_size);
          last_out_power = calculate_power(last_outputs,filter_size);
          
          if(last_out_power < vad_threshold || (!use_vad)){
            //updating filter G for next sample
            for (i = 0; i < number_of_microphones-1; i++){
              block_power = calculate_power(block_matrix[i],filter_size);
              
              //calculating the mu for this microphone
              if(mu0*block_power/last_out_power < mu_max){
                this_mu = mu0/last_out_power;
              }else{
                this_mu = mu0/block_power;
              }
              
              if(std::isnan(this_mu) || std::isinf(this_mu))
                this_mu = 0.0;
              
              //applying mu to filter update
              for(k = 0; k < filter_size; ++k){
                filter[i][k] += this_mu*out[j]*block_matrix[i][k];
                
                //important NaN check
                //if this is not done, the NaN is distributed to future updates
                if(std::isnan(filter[i][k]))
                  filter[i][k] = 0.0;
              }
              
              //registering mu of the first microphone for future reference
              if(write_mu && i == 0){
                avg_mu += this_mu;
              }
            }
          }else if(write_mu){
            avg_mu = last_avg_mu;
          }
        }
        
        if(write_mu){
          fprintf(mu_file,"%f\n",avg_mu/nframes);
          last_avg_mu = avg_mu;
        }
        
        //Outputing to ROS
        output_to_rosjack (out, nframes, output_type);
        
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

void gsc_handle_params(ros::NodeHandle *n){
    std::string node_name = ros::this_node::getName();
    std::cout << "GSC ROS parameters: " << std::endl;
    
    if ((*n).getParam(node_name+"/use_vad",use_vad)){
        ROS_INFO("Use VAD: %d",use_vad);
        if(use_vad){
            if ((*n).getParam(node_name+"/vad_threshold",vad_threshold)){
                ROS_INFO("VAD Threshold: %f",vad_threshold);
            }else{
                vad_threshold = 0.1;
                ROS_WARN("VAD Threshold argument not found in ROS param server, using default value (%f).",vad_threshold);
            }
        }
    }else{
        use_vad = false;
        ROS_WARN("Use VAD argument not found in ROS param server, using default value (%d).",use_vad);
    }
    
    if ((*n).getParam(node_name+"/mu0",mu0)){
        ROS_INFO("Initial Mu: %f",mu0);
    }else{
        mu0 = 0.0005;
        ROS_WARN("Initial Mu argument not found in ROS param server, using default value (%f).",mu0);
    }
    
    if ((*n).getParam(node_name+"/mu_max",mu_max)){
        ROS_INFO("Maximum Mu: %f",mu_max);
    }else{
        mu_max = 0.01;
        ROS_WARN("Maximum Mu argument not found in ROS param server, using default value (%f).",mu_max);
    }
    
    if ((*n).getParam(node_name+"/filter_size",filter_size)){
        ROS_INFO("Filter size: %d",filter_size);
    }else{
        filter_size = 128;
        ROS_WARN("Filter size argument not found in ROS param server, using default value (%d).",filter_size);
    }
    
    if ((*n).getParam(node_name+"/write_mu",write_mu)){
        ROS_INFO("Write Mu: %d",write_mu);
        if(write_mu){
          strcpy(mu_file_path,home_path);
          strcat(mu_file_path,(char *)"/mu_behavior.txt");
          ROS_INFO("Write Mu file path: %s",mu_file_path);
          mu_file = fopen(mu_file_path,"w");
        }
    }else{
        write_mu = false;
        ROS_WARN("Write Mu argument not found in ROS param server, using default value (%d).",write_mu);
    }
}


int main (int argc, char *argv[]) {
    int i,j;
    /* ROS initialization*/
    ros::init(argc, argv, client_name);
    ros::NodeHandle n;
    handle_params(&n);
    gsc_handle_params(&n);
    
    ros::Subscriber theta_subscriber = n.subscribe("theta", 1000, theta_roscallback);
    
    
    /* create JACK agent */
    if(rosjack_create (ROSJACK_READ, &n, "jackaudio", client_name, number_of_microphones, jack_callback)){
        ROS_ERROR("JACK agent could not be created.\n");
        ros::shutdown();
        exit(1);
    }
    
    std::cout << "Pre-allocating space for internal buffers." << std::endl;
    prepare_overlap_and_add_bymic(); //fft_win is assinged here
    
    x_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    x_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    
    x_forward = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(x_time), reinterpret_cast<fftw_complex*>(x_fft), FFTW_FORWARD, FFTW_MEASURE);
    y_inverse = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(y_fft), reinterpret_cast<fftw_complex*>(y_time), FFTW_BACKWARD, FFTW_MEASURE);
    
    block_matrix = (rosjack_data **) malloc (sizeof(rosjack_data*)*(number_of_microphones-1));
    filter = (rosjack_data **) malloc (sizeof(rosjack_data*)*(number_of_microphones-1));
    for (i = 0; i < number_of_microphones-1; i++){
        block_matrix[i] = (rosjack_data *) calloc (filter_size,sizeof(rosjack_data));
        filter[i] = (rosjack_data *) calloc (filter_size,sizeof(rosjack_data));
    }
    last_outputs = (rosjack_data *) calloc (filter_size,sizeof(rosjack_data));
    
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
