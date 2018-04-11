
# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# define project information
project(ElastoPlasticQPot)

# automatically parse the version number
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/ElastoPlasticQPot/ElastoPlasticQPot.h" _ElastoPlasticQPot_version_header)
string(REGEX MATCH "define[ \t]+ELASTOPLASTICQPOT_WORLD_VERSION[ \t]+([0-9]+)" _ElastoPlasticQPot_world_version_match "${_ElastoPlasticQPot_version_header}")
set(ELASTOPLASTICQPOT_WORLD_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+ELASTOPLASTICQPOT_MAJOR_VERSION[ \t]+([0-9]+)" _ElastoPlasticQPot_major_version_match "${_ElastoPlasticQPot_version_header}")
set(ELASTOPLASTICQPOT_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+ELASTOPLASTICQPOT_MINOR_VERSION[ \t]+([0-9]+)" _ElastoPlasticQPot_minor_version_match "${_ElastoPlasticQPot_version_header}")
set(ELASTOPLASTICQPOT_MINOR_VERSION "${CMAKE_MATCH_1}")
set(ELASTOPLASTICQPOT_VERSION_NUMBER ${ELASTOPLASTICQPOT_WORLD_VERSION}.${ELASTOPLASTICQPOT_MAJOR_VERSION}.${ELASTOPLASTICQPOT_MINOR_VERSION})

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# set c++ standard -> c++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# set paths where 'ElastoPlasticQPot' will be installed
# - root
set(ELASTOPLASTICQPOT_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
# - headers
set(ELASTOPLASTICQPOT_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
# - headers
set(INCLUDE_INSTALL_DIR
  "${CMAKE_INSTALL_INCLUDEDIR}/ElastoPlasticQPot"
  CACHE PATH "Path, relative to CMAKE_PREFIX_PATH, where the ElastoPlasticQPot header files are installed"
)
# - CMake configuration
set(CMAKEPACKAGE_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/ElastoPlasticQPot"
  CACHE PATH "Path where ElastoPlasticQPotConfig.cmake is installed"
)
# - pkg-config
set(PKGCONFIG_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/pkgconfig"
  CACHE PATH "Path where ElastoPlasticQPot.pc is installed"
)

# list all header files
set(ELASTOPLASTICQPOT_HEADERS
  ElastoPlasticQPot/Cartesian2d.cpp
  ElastoPlasticQPot/Cartesian2d.h
  ElastoPlasticQPot/Cartesian2d_Cusp.cpp
  ElastoPlasticQPot/Cartesian2d_Elastic.cpp
  ElastoPlasticQPot/Cartesian2d_Matrix.cpp
  ElastoPlasticQPot/Cartesian2d_Smooth.cpp
  ElastoPlasticQPot/ElastoPlasticQPot.h
)

# Disable pkg-config for native Windows builds
if(NOT WIN32 OR NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(ELASTOPLASTICQPOT_BUILD_PKGCONFIG "Build pkg-config .pc file for ElastoPlasticQPot" ON)
endif()

# build pkg-config file -> fill in relevant values
if(ELASTOPLASTICQPOT_BUILD_PKGCONFIG)
  configure_file(ElastoPlasticQPot.pc.in ElastoPlasticQPot.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ElastoPlasticQPot.pc DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# configure CMake file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ElastoPlasticQPotConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/ElastoPlasticQPotConfig.cmake
  PATH_VARS ELASTOPLASTICQPOT_INCLUDE_DIR ELASTOPLASTICQPOT_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# - CMake-file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/ElastoPlasticQPotConfig.cmake
  DESTINATION
  ${CMAKEPACKAGE_INSTALL_DIR}
)
# - headers
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/${ELASTOPLASTICQPOT_HEADERS}
  DESTINATION
  ${INCLUDE_INSTALL_DIR}/ElastoPlasticQPot
)

# print information to screen
message(STATUS "")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in '${CMAKE_INSTALL_PREFIX}'")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake . -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
