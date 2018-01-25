#!/usr/bin/env python

import rospy
from time import sleep
from std_msgs.msg import Float32

mu = 0.01
past_SIR = -100.0
past_theta = 1.0
def SIRcallback(data):
    global pub
    global mu
    global past_SIR
    global past_theta
    
    rospy.loginfo(rospy.get_caller_id() + "I heard SIR %f", data.data)
    SIR = data.data
    theta = past_theta - mu*(SIR-past_SIR)
    rospy.loginfo("\t past SIR %f", past_SIR)
    rospy.loginfo("\t past theta %f", past_theta)
    rospy.loginfo("\t SIR %f", SIR)
    rospy.loginfo("\t publishing theta %f", theta)
    past_SIR = SIR
    past_theta = theta
    sleep(1)
    pub.publish(theta)
    
    

rospy.init_node('SIR2theta', anonymous=True)

try:
    pub = rospy.Publisher('theta', Float32, queue_size=10)
    rospy.Subscriber("SIR", Float32, SIRcallback)
except rospy.ROSInterruptException:
    pass

sleep(1)
pub.publish(past_theta)
rospy.spin()
