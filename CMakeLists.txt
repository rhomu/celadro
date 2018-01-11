# This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.0)
project(celadro)
set(CMAKE_CXX_STANDARD 14)

# default is release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# set warning level
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

################################################################################
# directories
################################################################################

# directories
include_directories(src)
add_subdirectory(src)

################################################################################
# target definitions
################################################################################

# use all files in src/
file(GLOB_RECURSE sources src/*.cpp src/*.hpp)

# main executable
add_executable(celadro ${sources})

################################################################################
# dependencies
################################################################################

# check for boost libraries
find_package(Boost 1.36.0 COMPONENTS program_options REQUIRED)

# cuda
find_package(CUDA)
if(CUDA_FOUND)
  file(GLOB_RECURSE cuda_sources src/*.cu src/*.h)
  list(APPEND CUDA_NVCC_FLAGS "-D_CUDA;-arch=sm_61;-std=c++11;-O2;-DVERBOSE")
  cuda_compile(cuda_objects ${cuda_sources})
  add_executable(celadro-cuda ${sources} ${cuda_objects})
  target_link_libraries(celadro-cuda PUBLIC ${Boost_LIBRARIES} ${CUDA_LIBRARIES})
  set_target_properties(celadro-cuda PROPERTIES APPEND COMPILE_FLAGS "-D_CUDA_ENABLED")
  target_include_directories(celadro-cuda PUBLIC "${CUDA_INCLUDE_DIRS}")
endif()

# openmp
find_package(OpenMP)
if(OPENMP_FOUND)
  set_target_properties(celadro PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  target_link_libraries(celadro ${Boost_LIBRARIES} ${OpenMP_CXX_FLAGS})
else()
  target_link_libraries(celadro PUBLIC ${Boost_LIBRARIES})
endif()

# add hydra support (boost is fucked up there)
option(HYDRA "Make linking work on hydra work (as of 2017)" OFF)
if(HYDRA)
  set_target_properties(celadro PROPERTIES LINK_FLAGS "${CMAKE_LINK_FLAGS} -L/usr/local/shared/boost/1.64.0-gcc5.4.0/lib")
endif()

################################################################################
# testing
################################################################################

################################################################################
# packaging
################################################################################
