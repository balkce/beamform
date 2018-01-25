/**
 * Functions to ease connection to JACK
 */

#include "rosjack.h"
void rosjack_handle_params(ros::NodeHandle *n){
	std::string node_name = ros::this_node::getName();
	std::cout << "ROS Node name: " << node_name << std::endl;
	
	if ((*n).getParam(node_name+"/output_type",output_type)){
		if (output_type >=0 && output_type < ROSJACK_OUT_ENUM){
			ROS_INFO("Output type: %s",ROSJACK_OUT_OUTPUT_TYPES[output_type]);
		}else{
			output_type = ROSJACK_OUT_BOTH;
			ROS_WARN("Invalid output type argument, using default value (%s).",ROSJACK_OUT_OUTPUT_TYPES[output_type]);
		}
	}else{
		output_type = ROSJACK_OUT_BOTH;
		ROS_WARN("Output type argument not found in ROS param server, using default value (%s).",ROSJACK_OUT_OUTPUT_TYPES[output_type]);
	}
	
	if ((*n).getParam(node_name+"/auto_connect",auto_connect)){
		ROS_INFO("Auto connect: %d",auto_connect);
	}else{
		auto_connect = true;
		ROS_WARN("Auto connect argument not found in ROS param server, using default value (%d).",auto_connect);
	}
}

void jack_shutdown (void *arg){
	exit (1);
}

int rosjack_create (int rosjack_type, ros::NodeHandle *n, const char *topic_name, const char *client_name, int input_number, int (*callback_function)(jack_nframes_t, void*)) {
	/* ROS stuff */
	signal(SIGINT, siginthandler);
	if(rosjack_type == ROSJACK_READ){
		rosjack_out = (*n).advertise<jack_msgs::JackAudio>(topic_name, 1000);
	}
	
	rosjack_handle_params(n);
	
	/* JACK initialization */
	int i;
	jack_num_inputs = input_number;
	printf ("Connecting to Jack Server...\n");
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	
	/* open a client connection to the JACK server */
	jack_client = jack_client_open (client_name, options, &status);
	if (jack_client == NULL){
		/* if connection failed, say why */
		printf ("jack_client_open() failed, status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			printf ("Unable to connect to JACK server.\n");
		}
		return 1;
	}
	
	/* if connection was successful, check if the name we proposed is not in use */
	if (status & JackNameNotUnique){
		client_name = jack_get_client_name(jack_client);
		printf ("Warning: other agent with our name is running, `%s' has been assigned to us.\n", client_name);
	}
	
	/* tell the JACK server to call 'jack_callback()' whenever there is work to be done. */
	jack_set_process_callback (jack_client, callback_function, 0);
	
	
	/* tell the JACK server to call 'jack_shutdown()' if it ever shuts down,
	   either entirely, or if it just decides to stop calling us. */
	jack_on_shutdown (jack_client, jack_shutdown, 0);
	
	
	/* display the current sample rate. */
	rosjack_window_size = jack_get_buffer_size (jack_client);
	printf ("JACK window size: %d\n", rosjack_window_size);
	rosjack_sample_rate = jack_get_sample_rate (jack_client);
	printf ("JACK sample rate: %d\n", rosjack_sample_rate);
	
	
	/* create the agent input ports */
	jack_input_port = (jack_port_t **)malloc(sizeof(jack_port_t *)*jack_num_inputs);
	char input_port_name[100];
	for(i = 0; i < jack_num_inputs; i++){
		sprintf(input_port_name,"input_%d",i+1);
		jack_input_port[i] = jack_port_register (jack_client, input_port_name, JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);

		/* check that the port were created succesfully */
		if ((jack_input_port[i] == NULL)) {
			printf("Could not create input port %s. Have we reached the maximum amount of JACK input ports?\n",input_port_name);
			return 1;
		}
	}

	jack_output_port = jack_port_register (jack_client, "output", JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
	if ((jack_output_port == NULL)) {
		printf("Could not create output port. Have we reached the maximum amount of JACK output ports?\n");
		return 1;
	}
	
	
	if(rosjack_type == ROSJACK_WRITE){
		ros2jack_buffer_size = jack_get_buffer_size (jack_client)*10;
		ros2jack_buffer = (rosjack_data *)malloc(sizeof(rosjack_data)*ros2jack_buffer_size);
	}

	/* Tell the JACK server that we are ready to roll.
	   Our jack_callback() callback will start running now. */
	if (jack_activate (jack_client)) {
		printf ("Cannot activate JACK agent.");
		return 1;
	}
	
	printf ("Agent activated.\n");
	
	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */
	const char **serverports_names;
	if(auto_connect){
		printf ("Connecting input ports... ");
		 
		/* Assign our ports to a server ports*/
		// Find possible output server port names
		serverports_names = jack_get_ports (jack_client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
		if (serverports_names == NULL) {
			printf("No available physical capture (server output) ports.\n");
			return 1;
		}
		// Connect the first available to our input port
		for(i = 0; i < jack_num_inputs; i++){
			if (jack_connect (jack_client, serverports_names[i], jack_port_name (jack_input_port[i]))) {
				printf("Cannot connect input port %s.\n",jack_port_name (jack_input_port[i]));
				ROS_WARN("Not connecting any more input ports, sticking with the ones that were connected.\n");
				break;
			}
		}
		// free serverports_names variable for reuse in next part of the code
		free (serverports_names);
	}
	
	// Find possible input server port names
	serverports_names = jack_get_ports (jack_client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		return 1;
	}
	// Connect the first available to our output port
	if (jack_connect (jack_client, jack_port_name (jack_output_port), serverports_names[0])) {
		printf("Cannot connect output port.\n");
		return 1;
	}
	// free serverports_names variable for reuse in next part of the code
	free (serverports_names);
	
	
	if(rosjack_type == ROSJACK_WRITE){
		rosjack_in = (*n).subscribe(topic_name, 1000, rosjack_roscallback);
	}

	printf ("done.\n");
	return 0;
}

void siginthandler(int sig){
	ROS_INFO("Closing JACK client.");
	jack_client_close (jack_client);
	
	ROS_INFO("Closing ROS node.");
	ros::shutdown();
}

void close_rosjack(){
	jack_client_close (jack_client);
}

void output_to_rosjack (rosjack_data *data, int data_length, int out_type){
	output_type = out_type;
	output_to_rosjack (data, data_length);
}
void output_to_rosjack (rosjack_data *data, int data_length){
	/* This may seem as too much code repetition, but it's quicker this way online */
	/* The enclosing brackets in the switch are necessary to avoid re-definition errors */
	
	switch(output_type){
		case ROSJACK_OUT_BOTH:{
				jack_msgs::JackAudio out;
				out.size = data_length;
 				ros::Time win_stamp((double)jack_last_frame_time(jack_client)/rosjack_sample_rate);
				out.header.stamp = win_stamp;
				rosjack_data *out_j = (rosjack_data *)jack_port_get_buffer (jack_output_port, data_length);
				
				for (int j = 0; j < data_length; j++){
					out.data.push_back(data[j]);
					out_j[j] = data[j];
				}
				rosjack_out.publish(out);
			}break;
		
		case ROSJACK_OUT_JACK:{
				rosjack_data *out_j = (rosjack_data *)jack_port_get_buffer (jack_output_port, data_length);
				
				for (int j = 0; j < data_length; j++){
					out_j[j] = data[j];
				}
			}break;
		
		case ROSJACK_OUT_ROS:{
				jack_msgs::JackAudio out;
				out.size = data_length;
 				ros::Time win_stamp((double)jack_last_frame_time(jack_client)/rosjack_sample_rate);
				out.header.stamp = win_stamp;
				for (int j = 0; j < data_length; j++){
					out.data.push_back(data[j]);
				}
				rosjack_out.publish(out);
			}break;
	}
}

rosjack_data ** input_from_rosjack (int data_length){
	int i;
	
	rosjack_data **data = (rosjack_data **)malloc(sizeof(rosjack_data *)*jack_num_inputs);
	for (i = 0; i < jack_num_inputs; i++){
		data[i] = (rosjack_data *)jack_port_get_buffer (jack_input_port[i], data_length);
	}
	
	return data;
}

void rosjack_roscallback(const jack_msgs::JackAudio::ConstPtr& msg){
	int msg_size = msg->size;
	
	jack_mtx.lock();
	for (int i = 0; i < msg_size; i++){
		ros2jack_buffer[ros2jack_buffer_size_w] = msg->data[i];
		ros2jack_buffer_size_w++;
		if(ros2jack_buffer_size_w > ros2jack_buffer_size)
			ros2jack_buffer_size_w = 0;
	}
	jack_mtx.unlock();
}

rosjack_data * input_from_ros2jack_buffer (int data_length){
	rosjack_data *out = (rosjack_data *)malloc(sizeof(rosjack_data)*data_length);
	
	jack_mtx.lock();
	for (int i = 0; i < data_length; i++){
		out[i] = ros2jack_buffer[ros2jack_buffer_size_r];
		ros2jack_buffer[ros2jack_buffer_size_r] = 0.0;
		
		ros2jack_buffer_size_r++;
		if(ros2jack_buffer_size_r > ros2jack_buffer_size)
			ros2jack_buffer_size_r = 0;
	}
	jack_mtx.unlock();
	
	return out;
}
