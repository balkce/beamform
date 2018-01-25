cmake_minimum_required(VERSION 2.8.3)
project(jack_msgs)

find_package(catkin REQUIRED COMPONENTS
  std_msgs
  message_generation
)

add_message_files(FILES JackAudio.msg)

generate_messages(DEPENDENCIES std_msgs)

catkin_package(CATKIN_DEPENDS message_runtime)

include_directories(${catkin_INCLUDE_DIRS})
