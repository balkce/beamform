/**
 * Online geometric source separation, used in:
 * J. Valin, S. Yamamoto, J. Rouat, F. Michaud, K. Nakadai and H. G. Okuno, "Robust Recognition of Simultaneous Speech by a Mobile Robot," in IEEE Transactions on Robotics, vol. 23, no. 4, pp. 742-752, Aug. 2007.
 * https://ieeexplore.ieee.org/document/4285864
 * 
 * And implemented in:
 * https://github.com/introlab/odas/blob/master/src/system/steer2demixing.c
 */

#include "rosjack.h"
#include <beamform/InterfTheta.h>
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
double **interf_delays;
double alpha;
std::vector<Eigen::MatrixXcd> weights;    //A in paper
std::vector<Eigen::MatrixXcd> weights_h;
std::vector<Eigen::MatrixXcd> sep_matrix; //W in paper

//GSS config
double freq_mag_threshold = 1.5;
double freq_max = 4000;
double freq_min = 400;
double out_amp = 4.5;
double mu = 0.01;
double lambda = 0.0;
double interf_angle_threshold = 5.0;

//reused buffers
Eigen::MatrixXcd in_fft(1,1);
Eigen::MatrixXcd E(1,1);
Eigen::MatrixXcd dj1(1,1);
Eigen::MatrixXcd dj2(1,1);
Eigen::MatrixXcd identity(1,1);
Eigen::MatrixXcd this_yf(1,1);

void update_weights(bool ini=false){
    calculate_delays(delays);
    
    int i,j,k;
    
    //weights for DOI
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            if(ini){
                for(j = 0; j < fft_win; j++){
                    weights[j](i,0) = 1.0; 
                }
            }
        }else{
            for(j = 0; j < fft_win; j++){
                weights[j](i,0) = std::exp(-M_I*(double)2*PI*freqs[j]*delays[i]); 
            }
        }
    }
    
    //weights for interferences
    for(k = 0; k < interference_angles.size(); k++){
        calculate_interf_delays(interf_delays[k],k,interference_angles[k]);
        for(i = 0; i < array_geometry.size(); i++){
            if (i == 0){
                if(ini){
                    for(j = 0; j < fft_win; j++){
                        weights[j](i,k+1) = 1.0; 
                    }
                }
            }else{
                for(j = 0; j < fft_win; j++){
                    weights[j](i,k+1) = std::exp(-M_I*(double)2*PI*freqs[j]*interf_delays[k][i]); 
                }
            }
        }
    }
    
    //other weights 
    for(j = 0; j < fft_win; j++){
        weights_h[j] = weights[j].adjoint();
        sep_matrix[j] = weights[j].adjoint();
    }
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
                //obtain y applying sep_matrix to x
                this_yf = sep_matrix[j]*in_fft.col(j);
                y_fft[j] = this_yf(0); //outputting only the first source of interest (changing this to all sources will require re-writing the whole ROS package)
                
                //calculate E factor (covariance matrix with diagonal removed) used to calculate dj1 gradient
                E = this_yf*this_yf.adjoint();
                E.diagonal() -= E.diagonal();
                
                //calculate alpha factor (magnitude normalization) used to calculate dj1 gradient
                alpha = in_fft.col(j).squaredNorm();
                alpha *= alpha;
                
                //calculate gradients for both decorrelation constraint (dj1) and geometric constraint (dj2)
                dj1 = 4 * (interference_angles.size()+1)     * (1/alpha) * (E * this_yf) * in_fft.col(j).adjoint();
                dj2 = 2 * (1/(interference_angles.size()+1)) * ((sep_matrix[j] * weights[j]) - identity) * weights_h[j];
                
                //update sep_matrix
                sep_matrix[j] = ((1-lambda*mu)*sep_matrix[j]) - mu * (dj1 + dj2); //lambda produces a slow reduction of the overall energy, setting to zero for the time being
                
            }else{
                y_fft[j] = in_fft(0,j)*0.01;
                this_yf.setZero();
                this_yf(0,0) = y_fft[j];
            }
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
        //applying wola to avoid discontinuities in the time domain
        out[j] *= hann_win_wola[j];
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

void gss_handle_params(ros::NodeHandle *n){
    std::string node_name = ros::this_node::getName();
    std::cout << "GSS ROS parameters: " << std::endl;
    
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
    
    if ((*n).getParam(node_name+"/mu",mu)){
        ROS_INFO("Adaptation rate (mu): %f",mu);
    }else{
        mu = 0.01;
        ROS_WARN("Adaptation rate (mu) argument not found in ROS param server, using default value (%f).",mu);
    }
    
    if ((*n).getParam(node_name+"/lambda",lambda)){
        ROS_INFO("Adaptation rate 2 (lambda): %f",lambda);
    }else{
        lambda = 0.0;
        ROS_WARN("Adaptation rate 2 (lambda) argument not found in ROS param server, using default value (%f).",lambda);
    }
    
    if ((*n).getParam(node_name+"/interf_angle_threshold",interf_angle_threshold)){
        ROS_INFO("Angular Threshold Between Interferences: %f",interf_angle_threshold);
    }else{
        interf_angle_threshold = 5.0;
        ROS_WARN("Angular Threshold Between Interferences argument not found in ROS param server, using default value (%f).",interf_angle_threshold);
    }
    
}

void free_interf_buffers(){
    int i;
    
    for(i = 0; i<interference_angles.size();i++){
        free(interf_delays[i]);
    }
    free(interf_delays);
    
    std::vector<Eigen::MatrixXcd>().swap(weights);
    std::vector<Eigen::MatrixXcd>().swap(weights_h);
    std::vector<Eigen::MatrixXcd>().swap(sep_matrix);
}

void allocate_interf_buffers(){
    //buffers to update: interf_delays, weights, sep_matrix (the last two could be recreated using update_weights), this_yf
    int i;
    
    interf_delays = (double **) malloc (sizeof(double*)*interference_angles.size());
    for(i = 0; i<interference_angles.size();i++){
        interf_delays[i] = (double *) malloc (sizeof(double)*number_of_microphones);
    }
    
    
    for(i = 0; i<fft_win;i++){
        Eigen::MatrixXcd this_weights(number_of_microphones,interference_angles.size()+1);
        Eigen::MatrixXcd this_weights_h(interference_angles.size()+1,number_of_microphones);
        Eigen::MatrixXcd this_sep_matrix(interference_angles.size()+1,number_of_microphones);
        
        this_weights.setZero();
        this_weights_h.setZero();
        this_sep_matrix.setZero();
        
        weights.push_back(this_weights);
        weights_h.push_back(this_weights);
        sep_matrix.push_back(this_sep_matrix);
    }
    
    this_yf.resize(interference_angles.size()+1,1);
    
    identity.resize(interference_angles.size()+1,interference_angles.size()+1);
    identity.setIdentity();
    
    dj1.resize(interference_angles.size()+1,number_of_microphones);
    dj2.resize(interference_angles.size()+1,number_of_microphones);
}

void interf_theta_roscallback(const beamform::InterfTheta::ConstPtr& msg){
    if (msg->id >= 1 && msg->id <= interference_angles.size()){
        ROS_INFO("Updating weights for interference %d to angle %f", msg->id, msg->angle);
        
        interference_angles[msg->id-1] = msg->angle;
        
        //check if the new angle is too close to another interference
        int i;
        for(i = 0; i<interference_angles.size();i++){
            if (i != (msg->id-1) && abs(interference_angles[i]-msg->angle) < interf_angle_threshold){
                ROS_INFO("New angle of interference %d is too close to the angle of interference %d: %f", msg->id, i+1, interference_angles[i]);
                ROS_INFO("Removing interference %d ", msg->id);
                
                READY = false;
                millisleep(30);
                free_interf_buffers();
                interference_angles.erase(interference_angles.begin()+msg->id-1);
                allocate_interf_buffers();
                READY = true;
                break;
            }
        }
        
        update_weights();
    }else if(msg->id > interference_angles.size()){
        ROS_INFO("Interference id (%d) is greater than current number of interferences (%d)", msg->id, (int)interference_angles.size());
        ROS_INFO("Assuming that it is a new interfernce...");
        
        int i;
        for(i = 0; i<interference_angles.size();i++){
            if (abs(interference_angles[i]-msg->angle) < interf_angle_threshold){
                ROS_INFO("Angle of new interference %d is too close to the angle of interference %d: %f", ((int)interference_angles.size())+1, i+1, interference_angles[i]);
                ROS_INFO("Not adding interference %d ", msg->id);
                break;
            }
        }
        
        if (i == interference_angles.size()){
            ROS_INFO("Adding new interference %d in angle %f", ((int)interference_angles.size())+1, msg->angle);
            READY = false;
            millisleep(30);
            free_interf_buffers();
            interference_angles.push_back(msg->angle);
            allocate_interf_buffers();
            READY = true;
            
            update_weights();
        }
    }else{
        ROS_INFO("Invalid interference id (%d), should be >= 1", msg->id);
    }
}

int main (int argc, char *argv[]) {
    int i;
    /* ROS initialization*/
    ros::init(argc, argv, client_name);
    ros::NodeHandle n;
    handle_params(&n);
    gss_handle_params(&n);
    
    ros::Subscriber theta_subscriber = n.subscribe("theta", 1000, theta_roscallback);
    ros::Subscriber interf_theta_subscriber = n.subscribe("theta_interference", 1000, interf_theta_roscallback);
    
    
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
    
    in_fft.resize(number_of_microphones,fft_win);
    
    E.resize(number_of_microphones,number_of_microphones);
    
    allocate_interf_buffers();
    
    update_weights(true);
    
    READY = true;
    
    ROS_INFO("Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    ros::spin();
    
    /* apparently ROS requires the use of exit, instead of return for final cleaning up */
    exit(0);
}
