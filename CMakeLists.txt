cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Debug)


# --- some definitions ---
if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)


# --- project ---
project(tethex)


# --- useful variables ---
set(TESTING 1 CACHE BOOL "Do you want to launch testing procedure before using the program?")
set(DELETE_SIMPLICES 0 CACHE BOOL "Delete simplices after converting them into bricks (to economize some memory)")


# --- write the variables initialized here to config.h file ---
message("Options: TESTING = ${TESTING}, DELETE_SIMPLICES = ${DELETE_SIMPLICES}")
configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h")
include_directories(${PROJECT_BINARY_DIR}) # because config.h file is in binary directory we need to include it


# --- sources and headers ---
set(SRC_LIST main.cpp tethex.cpp)
set(HDR_LIST tethex.h)


# --- Google Test ---
if(TESTING MATCHES 1)
  message("ATTENTION! This is a version with testing features!")
  set(GTEST_DIR /home/artemiev/projects/gtest-1.6.0 CACHE PATH "Path to Google Test top level directory")
  find_library(GTEST_LIB gtest ${GTEST_DIR}/lib/.libs)
  if(GTEST_LIB)
    message("GTest library is ${GTEST_LIB}")
  else(GTEST_LIB)
    message(FATAL_ERROR "GTest was not found!")
  endif(GTEST_LIB)
  include_directories(${GTEST_DIR}/include)
endif(TESTING MATCHES 1)


# --- make ---
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
if(TESTING MATCHES 1)
  target_link_libraries(${PROJECT_NAME} ${GTEST_LIB} pthread)
endif(TESTING MATCHES 1)

