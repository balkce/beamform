#!/usr/bin/env python

import rospy
from time import sleep
from std_msgs.msg import Float32
from jack_msgs.msg import JackAudio
import math as m
import numpy as np
from collections import deque

num_win = 50
vad_threshold = 0.001 #0.0001
hist_bins_range = 8
mu = 25

num_win_i = 0
past_energy = -100.0
past_theta = 0.0
windows = deque([]);

hist_bins = []

def get_energy_from_list(data_list):
    data_list_squared = [i ** 2 for i in data_list]
    energy = m.sqrt(sum(data_list_squared)/len(data_list_squared))
    
    return energy

def get_energy_from_deque(data_deque):
    global hist_bins
    data_list = [item for sublist in list(data_deque) for item in sublist]
    
    data_np = np.abs(np.array(data_list))
    
    if(len(hist_bins) > 0):
        data_hist_values,data_hist_bins = np.histogram(data_np,hist_bins)
    else:
        data_hist_values,data_hist_bins = np.histogram(data_np,'fd') # Freedman Diaconis Estimator of histogram bins: Robust (resilient to outliers) estimator that takes into account data variability and data size.
        hist_bins = data_hist_bins
    
    ## mu : 25
    #max_hist_bin_i = np.argmax(data_hist_values)
    #hist_bins_range = int(float(len(data_hist_values)/20))
    #hist_bins_to_use = np.arange(np.maximum(0,max_hist_bin_i-hist_bins_range),np.minimum(len(data_hist_values),max_hist_bin_i+hist_bins_range+1))
    #energy = float(sum(data_hist_values[hist_bins_to_use]))/len(data_list)
    
    ## mu : 25
    #middle_hist_bin_i = int(float(len(data_hist_values)/10))
    #hist_bins_range = int(float(len(data_hist_values)/10))
    #hist_bins_to_use = np.arange(np.maximum(0,middle_hist_bin_i-hist_bins_range),np.minimum(len(data_hist_values),middle_hist_bin_i+hist_bins_range+1))
    #energy = float(sum(data_hist_values[hist_bins_to_use]))/len(data_list)
    
    ## mu : 25
    data_hist_values_p = data_hist_values.astype(float)/len(data_list)
    energy = np.sum(data_hist_bins[0:-1]*data_hist_values_p) #expected value
    
    ## mu : 25
    #energy = m.sqrt(np.mean(data_np ** 2))
    
    return energy

def energycallback(data):
    global pub
    global mu
    global past_energy
    global past_theta
    global num_win
    global num_win_i
    global windows
    
    this_win = list(data.data)
    this_win_energy = get_energy_from_list(this_win)
    if (this_win_energy >= vad_threshold):
        if(num_win_i < num_win):
            windows.append(this_win)
            num_win_i += 1
            
        else:
            windows.popleft()
            windows.append(this_win)
            
            if(past_energy == -100.0):
                past_energy = get_energy_from_deque(windows)
            
            energy = get_energy_from_deque(windows)
            
            theta = past_theta + mu*(energy-past_energy) # gradient ascent (the plus sign is important)
            if(theta > 180):
                theta = theta-360
            elif(theta < -180):
                theta = theta+360
            
            rospy.loginfo("---")
            rospy.loginfo("\t past energy %f", past_energy)
            rospy.loginfo("\t      energy %f", energy)
            rospy.loginfo("\t        diff %f", energy-past_energy)
            rospy.loginfo("\t        step %f", mu*(energy-past_energy))
            rospy.loginfo("\t    past theta %f", past_theta)
            rospy.loginfo("\t current theta %f", theta)
            pub.publish(theta)
            
            past_energy = energy
            past_theta = theta

rospy.loginfo("Starting energy2theta node...")

rospy.init_node('energy2theta', anonymous=True)
initial_angle_param_name = '/beamform/initial_angle'
if rospy.has_param(initial_angle_param_name):
    past_theta = rospy.get_param(initial_angle_param_name)
    rospy.loginfo("\t Initial angle: %f", past_theta)
else:
    rospy.error("Parameter /beamform/initial_angle not found. Is beamform node running?")
    exit()

try:
    pub = rospy.Publisher('theta', Float32, queue_size=10)
    rospy.Subscriber("jackaudio", JackAudio, energycallback)
except rospy.ROSInterruptException:
    pass

rospy.spin()
