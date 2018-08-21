cmake_minimum_required(VERSION 2.8.3)
project(beamform)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wno-narrowing") #necessary for the use of mutex

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  jack_msgs
)
find_package(cmake_modules REQUIRED) # for Eigen
find_package(Eigen REQUIRED)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES beamform
 CATKIN_DEPENDS roscpp rospy std_msgs jack_msgs
#  DEPENDS system_lib
)

INCLUDE(FindPkgConfig)
pkg_check_modules (GLIB2 glib-2.0)
pkg_check_modules (JACK jack)
pkg_check_modules (FFTW fftw3)
pkg_check_modules (EIGEN eigen)
pkg_check_modules (SNDFILE sndfile)

include_directories(
  src/rosjack
  ${GLIB2_INCLUDE_DIRS}
  ${JACK_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  ${SNDFILE_INCLUDE_DIRS}
)

add_library(rosjack
  src/rosjack/rosjack.cpp
)

add_executable(das src/das.cpp)
add_executable(das_time src/das_time.cpp)
add_executable(mvdr src/mvdr.cpp)
add_executable(beamform src/beamform.cpp)
add_executable(rosjack_write src/jack_write.cpp)
add_executable(rosjack_read src/jack_read.cpp)
add_executable(rosjack_ref src/jack_ref.cpp)
add_executable(phase src/phase.cpp)
add_executable(gsc src/gsc.cpp)

target_link_libraries(rosjack
  ${JACK_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  ${catkin_LIBRARIES}
)
add_dependencies(rosjack ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

target_link_libraries(das
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(das_time
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(mvdr
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(beamform
  rosjack
  ${catkin_LIBRARIES}
)

target_link_libraries(rosjack_read
  rosjack
  ${catkin_LIBRARIES}
)

target_link_libraries(rosjack_write
  rosjack
  ${catkin_LIBRARIES}
)

target_link_libraries(rosjack_ref
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(phase
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(gsc
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

