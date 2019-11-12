cmake_minimum_required(VERSION 2.8.3)
project(beamform)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wno-narrowing") #necessary for the use of mutex

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  jack_msgs
  message_generation
)
find_package(cmake_modules REQUIRED) # for Eigen
find_package(Eigen REQUIRED)

add_message_files(FILES InterfTheta.msg)

generate_messages(DEPENDENCIES std_msgs)

catkin_package(
  CATKIN_DEPENDS roscpp rospy std_msgs jack_msgs message_runtime
)

INCLUDE(FindPkgConfig)
pkg_check_modules (GLIB2 glib-2.0)
message( "Found glib-2.0: " ${GLIB2_FOUND} )
pkg_check_modules (JACK jack)
message( "Found jack: " ${JACK_FOUND} )
pkg_check_modules (FFTW fftw3)
message( "Found fftw3: " ${FFTW_FOUND} )
pkg_check_modules (EIGEN eigen)
message( "Found eigen: " ${EIGEN_FOUND} )
pkg_check_modules (SNDFILE sndfile)
message( "Found sndfile: " ${SNDFILE_FOUND} )
pkg_check_modules (SAMPLERATE samplerate)
message( "Found samplerate: " ${SAMPLERATE_FOUND} )

include_directories(
  src/rosjack
  ${GLIB2_INCLUDE_DIRS}
  ${JACK_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIRS}
  ${FFTW_INCLUDE_DIRS}
  ${SNDFILE_INCLUDE_DIRS}
  ${SAMPLERATE_INCLUDE_DIRS}
)

add_library(rosjack
  src/rosjack/rosjack.cpp
)

add_executable(das src/das.cpp)
add_executable(mvdr src/mvdr.cpp)
add_executable(rosjack_write src/jack_write.cpp)
add_executable(rosjack_read src/jack_read.cpp)
add_executable(rosjack_ref src/jack_ref.cpp)
add_executable(phase src/phase.cpp)
add_executable(phasempf src/phasempf.cpp)
add_executable(gsc src/gsc.cpp)
add_executable(lcmv src/lcmv.cpp)
add_executable(gss src/gss.cpp)
add_executable(mcra src/mcra.cpp)

target_link_libraries(rosjack
  ${JACK_LIBRARIES}
  ${SNDFILE_LIBRARIES}
  ${SAMPLERATE_LIBRARIES}
  ${catkin_LIBRARIES}
)
add_dependencies(rosjack ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

target_link_libraries(das
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(mvdr
  rosjack
  ${FFTW_LIBRARIES}
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

target_link_libraries(phasempf
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(gsc
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(lcmv
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(gss
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_link_libraries(mcra
  rosjack
  ${FFTW_LIBRARIES}
  ${catkin_LIBRARIES}
)


