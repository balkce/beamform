/**
 * Functions to ease connection to JACK
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <signal.h>
#include <sys/types.h>
#include <pwd.h>
#include <mutex>

#include <jack/jack.h>

#include <sndfile.h>

/*** ROS libraries ***/
#include "ros/ros.h"
#include "std_msgs/Float32.h"
#include "std_msgs/Header.h"
#include <jack_msgs/JackAudio.h>
/*** End ROS libraries ***/

#define ROSJACK_OUT_BOTH 0
#define ROSJACK_OUT_JACK 1
#define ROSJACK_OUT_ROS 2
#define ROSJACK_OUT_ENUM 3

#define ROSJACK_READ 0
#define ROSJACK_WRITE 1

char *rosjack_home_path;

//sndfile stuff
char audio_file_path[FILENAME_MAX];
SNDFILE * audio_file;
SF_INFO audio_info;
float *write_file_buffer;
int write_file_count;

const char *ROSJACK_OUT_OUTPUT_TYPES[] = {
  "ROSJACK_OUT_BOTH",
  "ROSJACK_OUT_JACK",
  "ROSJACK_OUT_ROS"
}; 

typedef jack_default_audio_sample_t rosjack_data;

jack_port_t    **jack_input_port;
jack_port_t    *jack_output_port;
jack_client_t  *jack_client;
int jack_num_inputs = 1;
int output_type;
ros::Publisher rosjack_out;
ros::Subscriber rosjack_in;

std::mutex jack_mtx;

bool auto_connect = true;
bool write_file = false;
bool write_xrun = false;

unsigned int xruns_count = 0;

rosjack_data *ros2jack_buffer;
unsigned int  ros2jack_buffer_size;
unsigned int ros2jack_buffer_size_r = 0;
unsigned int ros2jack_buffer_size_w = 0;
unsigned int rosjack_window_size = 0;
unsigned int rosjack_sample_rate = 0;

void rosjack_handle_params(ros::NodeHandle *n);
void jack_shutdown (void *arg);
int jack_xrun (void *arg);
void rosjack_roscallback(const jack_msgs::JackAudio::ConstPtr& msg);
int rosjack_create (int rosjack_type, ros::NodeHandle *n, const char *topic_name, const char *client_name, int input_number, int (*callback_function)(jack_nframes_t, void*));
void close_rosjack();
void siginthandler(int sig);
void output_to_rosjack (rosjack_data *data, int data_length, int output_type);
void output_to_rosjack (rosjack_data *data, int data_length);
rosjack_data ** input_from_rosjack (int data_length);
rosjack_data * input_from_ros2jack_buffer (int data_length);
