# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(lf_deal)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#---Create a shared library with geneated dictionary
#add_library(lf_deal SHARED lf_deal.cpp)
#target_link_libraries(lf_deal ${ROOT_LIBRARIES})

#---Create  a main program using the library
link_directories(/home/zqh/gsl/install/lib)
add_executable(lf_deal lf_deal.cpp ${sources} ${headers})
target_link_libraries(lf_deal ${ROOT_LIBRARIES} -lgsl -lgslcblas)

# Copy all scripts to the build directory
set(pals_SCRIPTS
  input.file
  )

foreach(_script ${pals_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
