#!/usr/bin/env python
# license removed for brevity
import rospy
from std_msgs.msg import String
from jack_msgs.msg import JackAudio
import threading
import time

startflag = False

this_win = []

def callback(data):
    global this_win
    global startflag
    
    # read from jackaudio
    this_win = list(data.data)
    #rospy.loginfo(rospy.get_caller_id() + "I heard %f : length %d", this_win[0], len(this_win))
    
    if not startflag:
        print (str(startflag))
        startflag = True

def thread_function(name):
    global pub
    global this_win
    global startflag
    
    #publish to speakerid topic
    while not startflag:
        time.sleep(1)
    
    rate = rospy.Rate(10) # 10hz
    while True:
        hello_str = "hello world %s" % this_win[0]
        rospy.loginfo(hello_str)
        pub.publish(hello_str)
        rate.sleep()

rospy.loginfo("Starting speakeridrest node...")

rospy.init_node('speakeridrest', anonymous=True)

try:
    pub = rospy.Publisher('speakerid', String, queue_size=10)
    rospy.Subscriber("jackaudio", JackAudio, callback)
    
    x = threading.Thread(target=thread_function, args=(1,))
    x.start()
    
except rospy.ROSInterruptException:
    pass

rospy.spin()
