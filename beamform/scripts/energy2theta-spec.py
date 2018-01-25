#!/usr/bin/env python

## requires python-scipy, python-numpy

import rospy
import message_filters
from time import sleep
from std_msgs.msg import Float32
from jack_msgs.msg import JackAudio
import math as m
import numpy as np
from collections import deque
from scipy import signal
import matplotlib
matplotlib.use('GTKAgg') 
import matplotlib.pylab as plt

energy_calc_method = 'history'

fs = 48000
num_win = 100 #500
vad_threshold = 0.001 #0.0001
fft_threshold = 0.00001
hist_bins_range = 8
mu = 5000

num_win_i = 0
past_energy = -100.0
past_theta = 0.0
windows = deque([]);

plt_i = 0
plt.axis([0, 300, -0.2, 0.2])
plt.ion()

def get_energy_from_list(data_list):
    data_list_squared = [i ** 2 for i in data_list]
    energy = m.sqrt(sum(data_list_squared)/len(data_list_squared))
    
    return energy

def get_energy_from_deque():
    global fft_threshold
    global fs
    global windows
    global mu
    global energy_calc_method
    global plt_i
    
    data_list = [item for sublist in list(windows) for item in sublist]
    
    data_np = np.array(data_list)
    
    energy = 0
    
    if (energy_calc_method == 'spectrogram'):
        ## this method measures the energy of the spectrogram of all the
        ## stored audio data using only high energy time-frequency bins
        
        ## it requires a constant energy from the sound sources to work well
        ## with a gradient descent
        ## it's not adequate for dynamic sources such as speech
        
        mu = 5000
        
        ## obtaining spectrogram of last num_win windows
        f, t, spec_data = signal.spectrogram(data_np, fs,nperseg=1024, noverlap=512, scaling='spectrum')
        
        #plotting
        #plt.pcolormesh(t, f, spec_data)
        #plt.pause(0.0001)
        
        # energy of the thresholded differential spectogram
        spec_data_filt = spec_data[spec_data > fft_threshold]
        energy = m.sqrt(np.mean(spec_data_filt))
        
    elif (energy_calc_method == 'history'):
        ## this method normalizes the energy considering past energy
        ## values, to constant-ify the search space and work well with
        ## a gradient descent
        
        mu = 10
        alpha = 1000
        
        past_values = np.array([np.sqrt(np.mean(np.array(win) ** 2)) for win in list(windows)])
        
        delta = past_values[-1]-np.mean(past_values)
        
        energy = past_values[-1]/(delta*alpha)
        
        plt.scatter(plt_i,past_values[-1],c='b')
        plt.scatter(plt_i,delta,c='r')
        plt.scatter(plt_i,energy,c='g')
        plt.pause(0.0001)
        plt_i+=1
        
    else:
        energy = -100.0
    
    
    if(m.isnan(energy)):
        energy = -100.0 #invalid
    
    return energy

def energycallback(jackaudio,jackaudio_ref):
    global pub
    global mu
    global past_energy
    global past_theta
    global num_win
    global num_win_i
    global windows
    global vad_threshold
    
    this_win = (np.array(list(jackaudio_ref.data))-np.array(list(jackaudio.data))).tolist()
    
    if(num_win_i < num_win):
        windows.append(this_win)
        num_win_i += 1
        
    else:
        if(num_win_i == num_win):
            print('Buffers full. Optimization can begin.')
            num_win_i += 1
        
        windows.popleft()
        windows.append(this_win)
        
        this_win_energy = get_energy_from_list(this_win)
        if (this_win_energy >= vad_threshold):
            if(past_energy == -100.0):
                past_energy = get_energy_from_deque()
            
            energy = get_energy_from_deque()
            
            if (energy > -100.0):
                theta = past_theta + mu*(energy-past_energy) # gradient descent (the minus sign is important)
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
                
                #comment next line to test get_energy_from_deque function with constant input
                pub.publish(theta)
                
                past_energy = energy
                past_theta = theta

rospy.loginfo("Starting energy2theta-sec node...")

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
    audiofil_sub = message_filters.Subscriber("jackaudio", JackAudio)
    audioref_sub = message_filters.Subscriber("jackaudio_ref", JackAudio)
    
    ts = message_filters.TimeSynchronizer([audiofil_sub, audioref_sub], 10)
    ts.registerCallback(energycallback)
except rospy.ROSInterruptException:
    pass

rospy.spin()
