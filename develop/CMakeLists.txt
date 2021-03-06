
cmake_minimum_required(VERSION 2.8.12)

project(test)

option(WARNINGS "Show build warnings" ON)

# basic compiler options
# ----------------------

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# switch on warnings
if(WARNINGS)
  if(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  endif()
endif()

# load packages
# -------------

# add custom paths
if(NOT "$ENV{INCLUDE_PATH}" STREQUAL "")
  string(REPLACE ":" ";" INCLUDE_LIST "$ENV{INCLUDE_PATH}")
  include_directories(SYSTEM ${INCLUDE_LIST})
endif()

# compile
# -------

add_executable(${PROJECT_NAME}
  main.cpp
  Cartesian2d.cpp
)
