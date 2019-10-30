#!/usr/bin/env python
# license removed for brevity
import rospy
from std_msgs.msg import String
from jack_msgs.msg import JackAudio
import threading
import time

this_win = []

count = 0

thread_alive = False

def callback(data):
    global this_win
    global count
    global thread_alive
    
    # read from jackaudio
    this_win = list(data.data)
    rospy.loginfo(rospy.get_caller_id() + "I heard %f : length %d", this_win[0], len(this_win))
    
    count = count + 1
    
    if(count > 10 and not thread_alive):
        x = threading.Thread(target=thread_function, args=(1,))
        x.start()
        count = 0

def thread_function(name):
    global pub
    global this_win
    global thread_alive
    
    thread_alive = True
    #while True:
    hello_str = "hello world %s" % this_win[0]
    rospy.loginfo(hello_str)
    pub.publish(hello_str)
    thread_alive = False


rospy.loginfo("Starting speakeridrest node...")

rospy.init_node('speakeridrest', anonymous=True)

try:
    pub = rospy.Publisher('speakerid', String, queue_size=10)
    rospy.Subscriber("jackaudio", JackAudio, callback)
    
    
except rospy.ROSInterruptException:
    pass

rospy.spin()
