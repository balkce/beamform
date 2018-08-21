/**
 * ROS agent that reads from microphone and outputs to ROS topic
 */

#include "rosjack.h"
#include "util.h"

int max_mic_past = -1;

int jack_callback (jack_nframes_t nframes, void *arg){
  int i,j;
  int max_mic_i = -1;
  double max_mic_energy = 0;
  double this_mic_energy = 0;
  
  //Inputing from ROS
  rosjack_data **in = input_from_rosjack (nframes);
  
  for (i = 0; i < jack_num_inputs; i++){
    this_mic_energy = 0;
    for(j = 0;j < nframes; j++)
      this_mic_energy += abs((double)(in[i][j]*100));
    
    //printf("%d %f      ",i,this_mic_energy);fflush(stdout);
    if (this_mic_energy > max_mic_energy || max_mic_past == -1){
      max_mic_energy = this_mic_energy;
      max_mic_i = i;
      
      if(max_mic_past == -1)
        max_mic_past = max_mic_i;
    }
  }
  
  if (max_mic_i == -1)
    max_mic_i = max_mic_past;
  
  //Outputing to ROS
  output_to_rosjack (in[max_mic_i], nframes, output_type);
  
  if(max_mic_past != max_mic_i && verbose)
    printf("Microphone %d now active.\n",max_mic_i+1);
  
  max_mic_past = max_mic_i;
  return 0;
}

int main (int argc, char *argv[]) {
  const char *client_name = "rosjack_read";
  
  /* ROS initialization*/
  ros::init(argc, argv, client_name);
  ros::NodeHandle n;
  handle_params(&n);
  
  /* create JACK agent */
  if(rosjack_create (ROSJACK_READ, &n, "jackaudio", client_name, number_of_microphones, jack_callback)){
    ROS_ERROR("JACK agent could not be created.\n");
    ros::shutdown();
    exit(1);
  }
  
  ROS_INFO("Beamform ROS node started.");
  
  /* keep running until stopped by the user */
  ros::spin();
  
  /* apparently ROS requires the use of exit, instead of return for final cleaning up */
  exit(0);
}
