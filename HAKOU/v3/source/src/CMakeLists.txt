#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/json)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

project(KSM C CXX)
cmake_minimum_required(VERSION 3.6)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


find_package(Eigen REQUIRED)
find_package(GenGetOpt REQUIRED)

LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/libs)
SET(json_LIB libjson.so)
#SET(json_LIB libjson_linux-gcc-4.2.1_libmt.a)

SET(SRC_LIST KSM_CalculationZip.cpp cmdline.c jsoncpp.cpp)

ADD_EXECUTABLE(KSM_CalculationZip ${SRC_LIST})

TARGET_LINK_LIBRARIES(KSM_CalculationZip)

set(CMAKE_BUILD_TYPE Release)

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O3 -msse2")        ## Optimize
endif()
