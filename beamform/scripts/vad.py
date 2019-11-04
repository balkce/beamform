#!/usr/bin/env python

import rospy
from time import sleep
from jack_msgs.msg import JackAudio
import math as m
import numpy as np

state_silence = False
state_active = False

tchange = 0.015
tvad = 0.02
ehist_len = 8
windows_passed_threshold = 5

ehist = np.zeros(ehist_len)
ehist_i = 0
enoise = 0.0
windows_passed = 0

def energycallback(data):
    global state_silence
    global state_active
    global tchange
    global tvad
    global ehist_len
    global ehist
    global ehist_i
    global enoise
    global windows_passed_threshold
    global windows_passed
    
    state_active 
    this_win = np.array(list(data.data))
    this_win_energy = np.absolute(this_win).mean()
    
    #checking for activity
    if (not state_silence and this_win_energy > enoise+tvad):
        windows_passed = 0
        state_active = True
        rospy.loginfo("VAD: active")
    else:
        state_active = False
        windows_passed = windows_passed + 1
    
    #checking for state change
    
    energy_mean = np.absolute(ehist).mean()
    
    if(state_silence and this_win_energy > energy_mean+tchange):
        state_silence = False
        enoise = energy_mean
        
        ehist = np.ones(ehist_len)*energy_mean
    elif(not state_silence and (this_win_energy < energy_mean-tchange or windows_passed > windows_passed_threshold)):
        windows_passed = 0
        state_silence = True
        
        ehist = np.ones(ehist_len)*enoise
    else:
        ehist[ehist_i] = this_win_energy
        ehist_i = ehist_i + 1
        if(ehist_i >= ehist_len):
            ehist_i = 0
    
    #print("Silence: " + str(state_silence) + " Active: " + str(state_active) + " this_win_energy: " + str(this_win_energy) + " enoise: " + str(enoise) + " energy_mean: " + str(energy_mean))


rospy.loginfo("Starting vad node...")

rospy.init_node('vad', anonymous=True)

try:
    rospy.Subscriber("jackaudio", JackAudio, energycallback)
except rospy.ROSInterruptException:
    pass

rospy.spin()
