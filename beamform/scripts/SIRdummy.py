#!/usr/bin/env python

import rospy
from std_msgs.msg import Float32

def thetacallback(data):
    global pub
    rospy.loginfo(rospy.get_caller_id() + "I heard theta %f", data.data)
    theta = data.data
    SIR = -(theta*theta)
    rospy.loginfo("\t publishing SIR %f", SIR)
    pub.publish(SIR)
    
    

rospy.init_node('SIR2dummy', anonymous=True)

try:
    pub = rospy.Publisher('SIR', Float32, queue_size=10)
    rospy.Subscriber("theta", Float32, thetacallback)
except rospy.ROSInterruptException:
    pass

rospy.spin()
