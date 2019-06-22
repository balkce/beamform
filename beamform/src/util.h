/*
    Some utilities for ROS handling and quality-of-life coding.
*/

//for variable handling
#include <map>
#include <vector>
#include <cmath>

//to measure execution time
#include <chrono>
typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

//for nanosleep function
#include <time.h>

//constants
#define PI 3.141592653589793238462643383279502884
std::complex<double> M_I(0,1);
double v_sound = 343;
double rad2deg = 180.0/PI;
double deg2rad = PI/180.0;
char *home_path;
const char *client_name = "beamform";

//global variables
bool verbose;
double angle;
int number_of_microphones = 1;
std::vector< std::map<std::string,double> > array_geometry;
std::vector< double > interference_angles;

//overlap-and-add stuff
unsigned int num_fftwindows = 4; //must be a power of 2
double *hann_win;
rosjack_data **in_buff;
rosjack_data **out_buff;
rosjack_data ***out_buff_mic;
rosjack_data **out_buff_mic_arg;
unsigned int out_buff_mic_arg_size;
unsigned int past_out_windows;
unsigned int out_buff_ini_shift;
unsigned int out_buff_last_shift;
unsigned int fft_win;

void handle_params(ros::NodeHandle *n){
    if ((home_path = getenv("HOME")) == NULL) {
        home_path = getpwuid(getuid())->pw_dir;
    }
    
    std::string node_name = ros::this_node::getName();
    std::cout << "ROS Node name: " << node_name << std::endl;
    std::cout << "Beamform ROS parameters: " << std::endl;
    
    if ((*n).getParam(node_name+"/verbose",verbose)){
        ROS_INFO("Verbose: %d",verbose);
    }else{
        verbose = false;
        ROS_WARN("Verbosity argument not found in ROS param server, using default value (%d).",verbose);
    }
    
    if ((*n).getParam(node_name+"/initial_angle",angle)){
        ROS_INFO("Initial angle: %f",angle);
    }else{
        angle = 0.0;
        ROS_WARN("Initial angle argument not found in ROS param server, using default value (%f).",angle);
    }
    
    int mic_number = 0;
    std::stringstream ss;
    ss << node_name+"/mic";
    ss << mic_number;
    std::string mic_param = ss.str();
    std::map<std::string,double> mic_map;
    
    while((*n).getParam(mic_param,mic_map)){
        mic_map["dist"] = sqrt(mic_map["x"]*mic_map["x"] + mic_map["y"]*mic_map["y"]);
        mic_map["angle"] = atan2(mic_map["y"],mic_map["x"]) * rad2deg;
        array_geometry.push_back(mic_map);
        ROS_INFO("Mic %d: (%f,%f) (%f,%f)",mic_number,mic_map["x"],mic_map["y"],mic_map["dist"],mic_map["angle"]);
        mic_number++;
        ss.str("");
        ss << node_name+"/mic";
        ss << mic_number;
        mic_param = ss.str();
    }
    
    int interf_number = 1;
    ss.str("");
    ss << node_name+"/angle_interf";
    ss << interf_number;
    std::string interf_param = ss.str();
    double interf_angle;
    
    while((*n).getParam(interf_param,interf_angle)){
        if (abs(interf_angle) <= 180){
            interference_angles.push_back(interf_angle);
            ROS_INFO("Interference %d at: %f",interf_number,interf_angle);
            interf_number++;
            ss.str("");
            ss << node_name+"/angle_interf";
            ss << interf_number;
            interf_param = ss.str();
        }else{
            break;
        }
    }
    
    //insert here code to consider first microphone as reference, ie. with coords (0,0)
    for(int i = 1; i < array_geometry.size(); i++){
        array_geometry[i]["x"] -= array_geometry[0]["x"];
        array_geometry[i]["y"] -= array_geometry[0]["y"];
    }
    
    //Doing some sane checks
    number_of_microphones = (int)array_geometry.size();
    ROS_INFO("Number of microphones: %d",number_of_microphones);
    
    std::cout << "Coordinates of microphones, relative to Mic0: " << std::endl;
    for(int i = 0; i < array_geometry.size(); i++){
        std::cout << "Mic #" << array_geometry[i]["id"] << std::endl;
        std::cout << "\t x: " << array_geometry[i]["x"] << std::endl;
        std::cout << "\t y: " << array_geometry[i]["y"] << std::endl;
        std::cout << "\t d: " << array_geometry[i]["dist"] << std::endl;
        std::cout << "\t a: " << array_geometry[i]["angle"] << std::endl;
    }
    std::cout << std::endl;
}

void calculate_delays(double *delay_buffer){
    int i,j;
    
    double this_dist = 0.0;
    double this_angle = 0.0;
    
    printf("New delays:\n");
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            //assuming first microphone as reference
            delay_buffer[i] = 0.0;
            printf("\t %d -> %f\n",i,delay_buffer[i]);
        }else{
            this_dist = array_geometry[i]["dist"];
            this_angle = array_geometry[i]["angle"]-angle;
            if(this_angle>180){
                this_angle -= 360;
            }else if(this_angle<-180){
                this_angle += 360;
            }
            
            delay_buffer[i] = this_dist*cos(this_angle*deg2rad)/(-v_sound);
            printf("\t %d -> %f\n",i,delay_buffer[i]);
        }
    }
}

void calculate_interf_delays(double *delay_buffer, int interf_id, double interf_angle){
    int i,j;
    
    double this_dist = 0.0;
    double this_angle = 0.0;
    
    printf("New delays for interference %d in: %f\n",interf_id+1,interf_angle);
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            //assuming first microphone as reference
            delay_buffer[i] = 0.0;
            printf("\t %d -> %f\n",i,delay_buffer[i]);
        }else{
            this_dist = array_geometry[i]["dist"];
            this_angle = array_geometry[i]["angle"]-interf_angle;
            if(this_angle>180){
                this_angle -= 360;
            }else if(this_angle<-180){
                this_angle += 360;
            }
            
            delay_buffer[i] = this_dist*cos(this_angle*deg2rad)/(-v_sound);
            printf("\t %d -> %f\n",i,delay_buffer[i]);
        }
    }
}

void calculate_frequency_vector(double *freq_buffer, unsigned int freq_buffer_size){
    int i;
    
    freq_buffer[0] = 0.0;
    for(i = 0; i<(freq_buffer_size/2)-1;i++){
        freq_buffer[i+1] = ((double)(i+1)/(double)freq_buffer_size)*((double)rosjack_sample_rate);
        freq_buffer[freq_buffer_size-1-i] = -((double)(i+1)/(double)freq_buffer_size)*((double)rosjack_sample_rate);
    }
    freq_buffer[(freq_buffer_size/2)-1] = ((double)rosjack_sample_rate)/2;
}

double hann(unsigned int buffer_i, unsigned int buffer_size){
    return 0.5 - 0.5*cos(2*PI*buffer_i/(buffer_size-1));
}

double* create_hann_winn (unsigned int h_size){
    static double *h = (double *) malloc(sizeof(double) * h_size);
    for (int i = 0; i < h_size; i++){
        h[i] = hann(i, h_size);
    }
    return h;
}

//fft_win is assigned here
//run before allocating buffers any other buffers
void prepare_overlap_and_add(){
    int i;
    
    fft_win = rosjack_window_size*num_fftwindows;
    past_out_windows = (int)(num_fftwindows/2)+1;
    out_buff_ini_shift = (int)((double)fft_win*3/4)-(int)((double)rosjack_window_size/2);//rosjack_window_size*2.5;
    out_buff_last_shift = (int)((double)fft_win/4)-(int)((double)rosjack_window_size/2);
    
    printf("Overlap Info:\n");
    printf("\t Start index of first window: %d\n",out_buff_ini_shift);
    printf("\t Start index of last window : %d\n",out_buff_last_shift);
    
    hann_win = create_hann_winn (fft_win);
    
    in_buff = (rosjack_data **) malloc (sizeof(rosjack_data*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; i++){
        in_buff[i] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
    }
    
    out_buff = (rosjack_data **) malloc (sizeof(rosjack_data*)*past_out_windows);
    for (i = 0; i < past_out_windows; i++){
        out_buff[i] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
    }
}

void do_overlap(rosjack_data **in, rosjack_data *out, jack_nframes_t nframes, void (*weight_func)(rosjack_data **, rosjack_data *)){
    int i,j;
    
    for (i = 0; i < number_of_microphones; i++){
        //appending this window to input buffer
        for(j = 0; j < nframes; j++)
            in_buff[i][j+(int)(nframes*(num_fftwindows-1))] = in[i][j];
    }
    
    //applying weights and storing the filter output in out_buff
    (*weight_func)(in_buff,out_buff[past_out_windows-1]);
    
    //doing overlap and storing in output
    for(j = 0; j < nframes; j++)
        out[j] = (out_buff[0][j+out_buff_ini_shift] + out_buff[past_out_windows-1][j+out_buff_last_shift]);
    
    //shifting input buffer one window
    for (i = 0; i < number_of_microphones; i++){
        for(j = 0;j < fft_win-nframes; j++)
            in_buff[i][j] = in_buff[i][j+nframes];
    }
    
    //shifting output buffer one window
    rosjack_data *out_bff_tmp = out_buff[0];
    for (i = 0; i < past_out_windows-1; i++){
        out_buff[i] = out_buff[i+1];
    }
    out_buff[past_out_windows-1] = out_bff_tmp;
}

//fft_win is assinged here
//run before allocating buffers any other buffers
void prepare_overlap_and_add_bymic(){
    int i,j;
    
    fft_win = rosjack_window_size*num_fftwindows;
    past_out_windows = (int)(num_fftwindows/2)+1;
    out_buff_ini_shift = (int)((double)fft_win*3/4)-(int)((double)rosjack_window_size/2);//rosjack_window_size*2.5;
    out_buff_last_shift = (int)((double)fft_win/4)-(int)((double)rosjack_window_size/2);
    
    printf("Overlap By Mic, Info:\n");
    printf("\t Start index of first window: %d\n",out_buff_ini_shift);
    printf("\t Start index of last window : %d\n",out_buff_last_shift);
    
    hann_win = (double *) malloc(sizeof(double) * fft_win);
    for (i = 0; i < fft_win; i++){
        hann_win[i] = hann(i, fft_win);
    }
    
    in_buff = (rosjack_data **) malloc (sizeof(rosjack_data*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; i++){
        in_buff[i] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
    }
    
    out_buff_mic = (rosjack_data ***) malloc (sizeof(rosjack_data**)*number_of_microphones);
    for (i = 0; i < number_of_microphones; i++){
        out_buff_mic[i] = (rosjack_data **) malloc (sizeof(rosjack_data*)*past_out_windows);
        for (j = 0; j < past_out_windows; j++){
            out_buff_mic[i][j] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
        }
    }
}

void do_overlap_bymic(rosjack_data **in, rosjack_data **out, jack_nframes_t nframes, void (*weight_func)(rosjack_data *, rosjack_data *, int)){
    int i,j;
    
    //doing overlap and add
    for (i = 0; i < number_of_microphones; i++){
      //appending this window to input buffer
      for(j = 0; j < nframes; j++)
          in_buff[i][j+(int)(nframes*(num_fftwindows-1))] = in[i][j];
      
      //applying weights and storing the filter output in out_buff_mic[i]
      (*weight_func)(in_buff[i],out_buff_mic[i][past_out_windows-1],i);
      
      //doing overlap and storing in output
      for(j = 0; j < nframes; j++)
          out[i][j] = (out_buff_mic[i][0][j+out_buff_ini_shift] + out_buff_mic[i][past_out_windows-1][j+out_buff_last_shift]);
      
      //shifting output buffer one window
      rosjack_data *out_bff_tmp = out_buff_mic[i][0];
      for (j = 0; j < past_out_windows-1; j++){
          out_buff_mic[i][j] = out_buff_mic[i][j+1];
      }
      out_buff_mic[i][past_out_windows-1] = out_bff_tmp;
    }
    
    //shifting in_overlap buffer
    for (i = 0; i < number_of_microphones; i++){
        for(j = 0;j < fft_win-nframes; j++)
            in_buff[i][j] = in_buff[i][j+nframes];
    }
}
//fft_win is assinged here
//run before allocating buffers any other buffers
void prepare_overlap_and_add_multi(int out_channels){
    int i,j;
    
    out_buff_mic_arg_size = out_channels;
    
    fft_win = rosjack_window_size*num_fftwindows;
    past_out_windows = (int)(num_fftwindows/2)+1;
    out_buff_ini_shift = (int)((double)fft_win*3/4)-(int)((double)rosjack_window_size/2);//rosjack_window_size*2.5;
    out_buff_last_shift = (int)((double)fft_win/4)-(int)((double)rosjack_window_size/2);
    
    printf("Overlap Multi, Info:\n");
    printf("\t Start index of first window: %d\n",out_buff_ini_shift);
    printf("\t Start index of last window : %d\n",out_buff_last_shift);
    
    hann_win = (double *) malloc(sizeof(double) * fft_win);
    for (i = 0; i < fft_win; i++){
        hann_win[i] = hann(i, fft_win);
    }
    
    in_buff = (rosjack_data **) malloc (sizeof(rosjack_data*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; i++){
        in_buff[i] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
    }
    
    out_buff_mic = (rosjack_data ***) malloc (sizeof(rosjack_data**)*out_channels);
    out_buff_mic_arg = (rosjack_data **) malloc (sizeof(rosjack_data*)*out_channels);
    for (i = 0; i < out_channels; i++){
        out_buff_mic[i] = (rosjack_data **) malloc (sizeof(rosjack_data*)*past_out_windows);
        out_buff_mic_arg[i] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
        for (j = 0; j < past_out_windows; j++){
            out_buff_mic[i][j] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
        }
    }
}

void do_overlap_multi(rosjack_data **in, rosjack_data **out, jack_nframes_t nframes, void (*weight_func)(rosjack_data **, rosjack_data **)){
    int i,j;
    
    for (i = 0; i < number_of_microphones; i++){
        //appending this window to input buffer
        for(j = 0; j < nframes; j++)
            in_buff[i][j+(int)(nframes*(num_fftwindows-1))] = in[i][j];
    }
    
    //applying weights and storing the filter output in out_buff_mic_arg
    (*weight_func)(in_buff,out_buff_mic_arg);
    
    //copying out_buff_mic_arg out_buff_mic[:][past_out_windows-1]
    for (i = 0; i < out_buff_mic_arg_size; i++){
        for (j = 0; j < fft_win; j++){
            out_buff_mic[i][past_out_windows-1][j] = out_buff_mic_arg[i][j];
        }
    }
    
    //doing overlap and storing in output
    for (i = 0; i < out_buff_mic_arg_size; i++){
        for(j = 0; j < nframes; j++)
            out[i][j] = (out_buff_mic[i][0][j+out_buff_ini_shift] + out_buff_mic[i][past_out_windows-1][j+out_buff_last_shift]);
    }
    
    //shifting input buffer one window
    for (i = 0; i < number_of_microphones; i++){
        for(j = 0;j < fft_win-nframes; j++)
            in_buff[i][j] = in_buff[i][j+nframes];
    }
    
    //shifting output buffer one window
    for (i = 0; i < out_buff_mic_arg_size; i++){
        rosjack_data *out_bff_tmp = out_buff_mic[i][0];
        for (j = 0; j < past_out_windows-1; j++){
            out_buff_mic[i][j] = out_buff_mic[i][j+1];
        }
        out_buff_mic[i][past_out_windows-1] = out_bff_tmp;
    }
}

//Sleep function in milliseconds
void millisleep(int milli){
     struct timespec st = {0};
     st.tv_sec = (milli/1000);
     st.tv_nsec = (milli%1000)*1000000L;
     nanosleep(&st, NULL);
}
