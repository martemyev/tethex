cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Debug)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

set(SRC_LIST main.cpp tethex.cpp)
set(HDR_LIST tethex.h)

project(tethex)
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})

