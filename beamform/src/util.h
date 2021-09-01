/*
    Some utilities for ROS handling and quality-of-life coding.
*/

//for variable handling
#include <map>
#include <vector>
#include <cmath>

//for rinbbuffer handling
#include <jack/ringbuffer.h>

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
// it carries out a 50% hop across the board
double *hann_win;
double hann_compensate;
jack_ringbuffer_t **in_buff;
rosjack_data **out_buff;
rosjack_data ***out_buff_mic;
rosjack_data **out_buff_mic_arg;
rosjack_data *out_buff_switch_tmp;
unsigned int out_buff_mic_arg_size;
unsigned int fft_win;
jack_ringbuffer_data_t *readring;
unsigned int bytes_nframes;

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
    for(int i = 1; i < array_geometry.size(); ++i){
        array_geometry[i]["x"] -= array_geometry[0]["x"];
        array_geometry[i]["y"] -= array_geometry[0]["y"];
    }
    
    //Doing some sane checks
    number_of_microphones = (int)array_geometry.size();
    ROS_INFO("Number of microphones: %d",number_of_microphones);
    
    std::cout << "Coordinates of microphones, relative to Mic0: " << std::endl;
    for(int i = 0; i < array_geometry.size(); ++i){
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
    for(i = 0; i < array_geometry.size(); ++i){
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
    for(i = 0; i < array_geometry.size(); ++i){
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
    for(i = 0; i<(freq_buffer_size/2)-1;++i){
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
    for (int i = 0; i < h_size; ++i){
        h[i] = sqrt(hann(i, h_size));
    }
    return h;
}

void prepare_hann(){
    int i;
    hann_win = create_hann_winn (fft_win);
    
    double hann_sum = 0;
    for (i = 0; i < fft_win; ++i){
        hann_sum += hann_win[i]*hann_win[i];
    }
    hann_compensate = ((double)rosjack_window_size)/hann_sum;
    std::cout << "Hann compensation: " << hann_compensate << std::endl;
}

void overlap_and_add_prepare_input(jack_ringbuffer_t *data, std::complex<double> *x){
    rosjack_data * buf;
    int buf_len, j;
    int i = 0;
    
    jack_ringbuffer_get_read_vector	(data,readring);

    // reconstructing current window from ring buffer data:
    //    it is built by two jack_ringbuffer_data_t-s, each with:
    //              len: bytes that are occupied
    //              buf: data to be read
    //    each jack_ringbuffer_data_t stores the 1st and 2nd "chunks"
    //    of data from the ringbuffer

    // extracting data from first "chunk" of ringbuffer
    buf = (rosjack_data*)(readring[0].buf);
    buf_len = (readring[0].len)/sizeof(jack_default_audio_sample_t);
    for(j=0; j < buf_len; ++i,++j)
      x[i] = buf[j]*hann_win[i]; //applying hann
    
    // extracting data from second "chunk" of ringbuffer
    buf = (rosjack_data*)(readring[1].buf);
    buf_len = (readring[1].len)/sizeof(jack_default_audio_sample_t);
    for(j=0; j < buf_len; ++i,++j)
      x[i] = buf[j]*hann_win[i]; //applying hann
}

void overlap_and_add_prepare_output(std::complex<double> *y, rosjack_data *out){
    int j;
    
    for (j = 0; j<fft_win; j++){
        // fftw3 does an unnormalized ifft that requires this normalization
        out[j] = real(y[j])/(double)fft_win;
        //applying wola to avoid discontinuities in the time domain
        out[j] *= hann_win[j];
    }
}

//fft_win is assigned here
//run before allocating buffers any other buffers
void prepare_overlap_and_add(){
    int i,j;
    rosjack_data zero_rosjack = 0.0;
    
    fft_win = rosjack_window_size*2;
    
    prepare_hann();
    
    in_buff = (jack_ringbuffer_t **) malloc (sizeof(jack_ringbuffer_t*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; ++i){
        // preparing the ring buffer
        //   jack_ringbuffer_t requires at least one sample between read and write pointers
        //   this is solved by using one addittional sample of memory
        //   if not, it will not provide the whole set of samples when reading from it
        //   which introduces problems when FFTW3 requires the full set
        in_buff[i] = jack_ringbuffer_create((fft_win+1)*sizeof(rosjack_data));
        jack_ringbuffer_reset(in_buff[i]); //resetting both read and write pointers
        
        // advancing and filling the ring buffer with zeros for one nframe
        for(j = 0; j < rosjack_window_size; ++j)
          jack_ringbuffer_write(in_buff[i],(char*)&zero_rosjack,sizeof(rosjack_data));
    }
  
    // other bits relevant to ring buffer
    bytes_nframes = rosjack_window_size*sizeof(rosjack_data);
    readring = (jack_ringbuffer_data_t *)malloc(2 * sizeof(jack_ringbuffer_data_t));
    
    out_buff = (rosjack_data **) malloc (sizeof(rosjack_data*)*2);
    out_buff[0] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
    out_buff[1] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
}

void do_overlap(rosjack_data **in, rosjack_data *out, jack_nframes_t nframes, void (*weight_func)(jack_ringbuffer_t **, rosjack_data *)){
    int i,j;
    
    for (i = 0; i < number_of_microphones; ++i){
        //adding this window to input buffer
        jack_ringbuffer_write(in_buff[i],(char*)in[i],bytes_nframes);
    }
    
    //applying weights and storing the filter output in out_buff
    (*weight_func)(in_buff,out_buff[1]);
    
    //doing overlap and storing in output
    for(j = 0; j < nframes; ++j)
        out[j] = (out_buff[0][j+nframes] + out_buff[1][j])*hann_compensate;
    
    //shifting input buffer one window
    for (i = 0; i < number_of_microphones; ++i){
        // advance ring buffer read pointer nframes
        jack_ringbuffer_read_advance(in_buff[i],bytes_nframes);
    }
    
    //shifting output buffer windows
    out_buff_switch_tmp = out_buff[0];
    out_buff[0] = out_buff[1];
    out_buff[1] = out_buff_switch_tmp;
}

//fft_win is assinged here
//run before allocating buffers any other buffers
void prepare_overlap_and_add_bymic(){
    int i,j;
    rosjack_data zero_rosjack = 0.0;

    fft_win = rosjack_window_size*2;
    
    prepare_hann();
    
    in_buff = (jack_ringbuffer_t **) malloc (sizeof(jack_ringbuffer_t*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; ++i){
        // preparing the ring buffer
        //   jack_ringbuffer_t requires at least one sample between read and write pointers
        //   this is solved by using one addittional sample of memory
        //   if not, it will not provide the whole set of samples when reading from it
        //   which introduces problems when FFTW3 requires the full set
        in_buff[i] = jack_ringbuffer_create((fft_win+1)*sizeof(rosjack_data));
        jack_ringbuffer_reset(in_buff[i]); //resetting both read and write pointers
        
        // advancing and filling the ring buffer with zeros for one nframe
        for(j = 0; j < rosjack_window_size; ++j)
          jack_ringbuffer_write(in_buff[i],(char*)&zero_rosjack,sizeof(rosjack_data));
    }
  
    // other bits relevant to ring buffer
    bytes_nframes = rosjack_window_size*sizeof(rosjack_data);
    readring = (jack_ringbuffer_data_t *)malloc(2 * sizeof(jack_ringbuffer_data_t));
    
    out_buff_mic = (rosjack_data ***) malloc (sizeof(rosjack_data**)*number_of_microphones);
    for (i = 0; i < number_of_microphones; ++i){
        out_buff_mic[i] = (rosjack_data **) malloc (sizeof(rosjack_data*)*2);
        out_buff_mic[i][0] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
        out_buff_mic[i][1] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
    }
}

void do_overlap_bymic(rosjack_data **in, rosjack_data **out, jack_nframes_t nframes, void (*weight_func)(jack_ringbuffer_t *, rosjack_data *, int)){
    int i,j;
    
    //doing overlap and add
    for (i = 0; i < number_of_microphones; ++i){
        //adding this window to input buffer
        jack_ringbuffer_write(in_buff[i],(char*)in[i],bytes_nframes);
        
        //applying weights and storing the filter output in out_buff_mic[i]
        (*weight_func)(in_buff[i],out_buff_mic[i][1],i);
        
        //doing overlap and storing in output
        for(j = 0; j < nframes; ++j)
            out[i][j] = (out_buff_mic[i][0][j+nframes] + out_buff_mic[i][1][j])*hann_compensate;
        
        //shifting output buffer one window
        out_buff_switch_tmp = out_buff_mic[i][0];
        out_buff_mic[i][0] = out_buff_mic[i][1];
        out_buff_mic[i][1] = out_buff_switch_tmp;
    }
    
    //shifting in_overlap buffer
    for (i = 0; i < number_of_microphones; ++i){
        // advance ring buffer read pointer nframes
        jack_ringbuffer_read_advance(in_buff[i],bytes_nframes);
    }
}
//fft_win is assinged here
//run before allocating buffers any other buffers
void prepare_overlap_and_add_multi(int out_channels){
    int i,j;
    rosjack_data zero_rosjack = 0.0;

    out_buff_mic_arg_size = out_channels;
    
    fft_win = rosjack_window_size*2;
    
    prepare_hann();
    
    in_buff = (jack_ringbuffer_t **) malloc (sizeof(jack_ringbuffer_t*)*number_of_microphones);
    for (i = 0; i < number_of_microphones; ++i){
        // preparing the ring buffer
        //   jack_ringbuffer_t requires at least one sample between read and write pointers
        //   this is solved by using one addittional sample of memory
        //   if not, it will not provide the whole set of samples when reading from it
        //   which introduces problems when FFTW3 requires the full set
        in_buff[i] = jack_ringbuffer_create((fft_win+1)*sizeof(rosjack_data));
        jack_ringbuffer_reset(in_buff[i]); //resetting both read and write pointers
        
        // advancing and filling the ring buffer with zeros for one nframe
        for(j = 0; j < rosjack_window_size; ++j)
          jack_ringbuffer_write(in_buff[i],(char*)&zero_rosjack,sizeof(rosjack_data));
    }
  
    // other bits relevant to ring buffer
    bytes_nframes = rosjack_window_size*sizeof(rosjack_data);
    readring = (jack_ringbuffer_data_t *)malloc(2 * sizeof(jack_ringbuffer_data_t));
    
    out_buff_mic = (rosjack_data ***) malloc (sizeof(rosjack_data**)*out_channels);
    out_buff_mic_arg = (rosjack_data **) malloc (sizeof(rosjack_data*)*out_channels);
    for (i = 0; i < out_channels; ++i){
        out_buff_mic[i] = (rosjack_data **) malloc (sizeof(rosjack_data*)*2);
        out_buff_mic_arg[i] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
        out_buff_mic[i][0] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
        out_buff_mic[i][1] = (rosjack_data *) calloc (fft_win,sizeof(rosjack_data));
    }
}

void do_overlap_multi(rosjack_data **in, rosjack_data **out, jack_nframes_t nframes, void (*weight_func)(jack_ringbuffer_t **, rosjack_data **)){
    int i,j;
    
    for (i = 0; i < number_of_microphones; ++i){
        //adding this window to input buffer
        jack_ringbuffer_write(in_buff[i],(char*)in[i],bytes_nframes);
    }
    
    //applying weights and storing the filter output in out_buff_mic_arg
    (*weight_func)(in_buff,out_buff_mic_arg);
    
    //copying out_buff_mic_arg to out_buff_mic[:][1]
    for (i = 0; i < out_buff_mic_arg_size; ++i){
        for (j = 0; j < fft_win; ++j){
            out_buff_mic[i][1][j] = out_buff_mic_arg[i][j];
        }
    }
    
    //doing overlap and storing in output
    for (i = 0; i < out_buff_mic_arg_size; ++i){
        for(j = 0; j < nframes; ++j)
            out[i][j] = (out_buff_mic[i][0][j+nframes] + out_buff_mic[i][1][j])*hann_compensate;
    }
    
    //shifting input buffer one window
    for (i = 0; i < number_of_microphones; ++i){
        // advance ring buffer read pointer nframes
        jack_ringbuffer_read_advance(in_buff[i],bytes_nframes);
    }
    
    //shifting output buffer one window
    for (i = 0; i < out_buff_mic_arg_size; ++i){
        out_buff_switch_tmp = out_buff_mic[i][0];
        out_buff_mic[i][0] = out_buff_mic[i][1];
        out_buff_mic[i][1] = out_buff_switch_tmp;
    }
}

//Sleep function in milliseconds
void millisleep(int milli){
     struct timespec st = {0};
     st.tv_sec = (milli/1000);
     st.tv_nsec = (milli%1000)*1000000L;
     nanosleep(&st, NULL);
}
