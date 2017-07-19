cmake_minimum_required(VERSION 2.8)
project(rosPE_link)

find_package(catkin REQUIRED COMPONENTS
   roscpp
   rospy
   std_msgs
   message_generation
   geometry_msgs
)

include_directories(
    ${OpenCV_INCLUDE_DIRS}
)

catkin_package(
  CATKIN_DEPENDS message_runtime
)


add_message_files(
  controlstack.msg
  posestack.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
  geometry_msgs
)

