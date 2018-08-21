/**
 * ROS agent that reads from ROS topic and outputs to speaker
 */

#include "rosjack.h"

int jack_callback (jack_nframes_t nframes, void *arg){
  output_to_rosjack (input_from_ros2jack_buffer (nframes), nframes, ROSJACK_OUT_JACK);
  return 0;
}


int main (int argc, char *argv[]) {
  const char *client_name = "rosjack_write";
  
  /* ROS initialization*/
  ros::init(argc, argv, client_name);
  ros::NodeHandle n;
  
  /* create JACK agent */
  if(rosjack_create (ROSJACK_WRITE,&n, "jackaudio", client_name, 0, jack_callback)){
    ROS_ERROR("JACK agent could not be created.\n");
    ros::shutdown();
    exit(1);
  }
  
  output_type = 2;
  auto_connect = 1;
  
  ROS_INFO("JackWrite ROS node started.");
  
  /* keep running until stopped by the user */
  ros::spin();
  
  /* apparently ROS requires the use of exit, instead of return for final cleaning up */
  exit (0);
}
