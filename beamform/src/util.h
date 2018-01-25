/*
    Some utilities for ROS handling and quality-of-life coding.
*/

//to measure execution time
#include <chrono>
typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

//constants
#define PI 3.141592653589793238462643383279502884
std::complex<double> M_I(0,1);
double v_sound = 343;
double rad2deg = 180.0/PI;
double deg2rad = PI/180.0;

//global variables
bool verbose;
double angle;
int number_of_microphones = 1;
std::vector< std::map<std::string,double> > array_geometry;


void handle_params(ros::NodeHandle *n){
    std::string node_name = ros::this_node::getName();
    std::cout << "ROS Node name: " << node_name << std::endl;
    
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
