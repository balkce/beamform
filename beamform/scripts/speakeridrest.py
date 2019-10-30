#!/usr/bin/env python
# license removed for brevity
import rospy
from std_msgs.msg import String
from jack_msgs.msg import JackAudio

def callback(data):
    global pub
    
    # read from jackaudio
    this_win = list(data.data)
    rospy.loginfo(rospy.get_caller_id() + "I heard %f : length %d", this_win[0], len(this_win))
    
    #publish to speakerid topic
    hello_str = "hello world %s" % rospy.get_time()
    rospy.loginfo(hello_str)
    pub.publish(hello_str)
    

rospy.loginfo("Starting speakeridrest node...")

rospy.init_node('speakeridrest', anonymous=True)

try:
    pub = rospy.Publisher('speakerid', String, queue_size=10)
    rospy.Subscriber("jackaudio", JackAudio, callback)
except rospy.ROSInterruptException:
    pass

rospy.spin()
