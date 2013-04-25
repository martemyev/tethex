cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Debug)


# --- some definitions ---
if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)


# --- project ---
project(tethex)


# --- set version ---
set(TETHEX_VERSION_MAJOR 0)
set(TETHEX_VERSION_MINOR 0)
set(TETHEX_VERSION_PATCH 1)


# --- useful variables ---
#set(TESTING 1 CACHE BOOL "Do you want to launch testing procedure before using the program?")
#set(DELETE_SIMPLICES 0 CACHE BOOL "Delete simplices after converting them into bricks (to economize some memory)")
option(TESTING "Do you want to launch testing procedure before using the program?" OFF)
option(DELETE_SIMPLICES "Delete simplices after converting them into bricks (to economize some memory)" ON)


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
if(TESTING)
  message("ATTENTION! This is a version with testing features!")
  set(GTEST_DIR /home/artemiev/projects/gtest-1.6.0 CACHE PATH "Path to Google Test top level directory")
  find_library(GTEST_LIB gtest ${GTEST_DIR}/lib/.libs)
  if(GTEST_LIB)
    message("GTest library is ${GTEST_LIB}")
  else(GTEST_LIB)
    message(FATAL_ERROR "GTest was not found!")
  endif(GTEST_LIB)
  include_directories(${GTEST_DIR}/include)
endif(TESTING)


# --- make ---
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
if(TESTING)
  target_link_libraries(${PROJECT_NAME} ${GTEST_LIB} pthread)
endif(TESTING)


# --- packaging of the library ---
set(CPACK_PACKAGE_VERSION_MAJOR "${TETHEX_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${TETHEX_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${TETHEX_VERSION_PATCH}")
# -- make package_source --
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "./*.user;build;doc;.git;.gitignore")
include(CPack)


# --- autodocumentation using doxygen ---
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(
    "${PROJECT_SOURCE_DIR}/doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/doxyfile")
  add_custom_target(doc COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)
else(DOXYGEN_FOUND)
  message("WARNING: Doxygen is not found - Autodocumentation will not be created")
endif(DOXYGEN_FOUND)
