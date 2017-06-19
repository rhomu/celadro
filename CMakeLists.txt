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
add_executable (celadro ${sources})

################################################################################
# dependencies
################################################################################

# check for boost libraries
find_package(Boost 1.36.0 COMPONENTS program_options REQUIRED)

target_link_libraries(celadro PUBLIC
  ${Boost_LIBRARIES}
)

# openmp
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

################################################################################
# testing
################################################################################

################################################################################
# packaging
################################################################################
