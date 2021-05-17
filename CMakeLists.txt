cmake_minimum_required(VERSION 3.17)
project(mtt_icm)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(PythonInterp 3)
find_package(PythonLibs 3)

# Install SFMT lib first (https://github.com/heavywatal/sfmt-class)
# Using homebrew: brew install heavywatal/tap/sfmt-class
find_package(sfmt)

# Tweak to make python discovery work (in my case it was pointing to the wrong install)
# First line doesn't seem to be necessary in the end, keeping it for the record
# set(PYTHON_INCLUDE_DIRS "/usr/local/Cellar/python@3.8/3.8.7_2/Frameworks/Python.framework/Versions/3.8/include/python3.8")
# Second line seems necessary only for the main target
# set(PYTHON_LIBRARY "/usr/local/Cellar/python@3.8/3.8.7_2/Frameworks/Python.framework/Versions/3.8/lib/libpython3.8.dylib")

# Standard boost installation
# Download boost https://www.boost.org/users/download/ and unzip
# Inside boost dir:
#    ./bootstrap.sh
#    ./b2
#    ./b2 install
# Other (brew / macports) install will have all chances to fail (at the time of writing)


# Components should be found in path with name libboost_NAME.dylib, eg: /usr/local/lib/libboost_numpy38.dylib
# This is additional to standard boost installation
# For both python38 and numpy38 to be produced, you have to:
# - be in the right python environment (conda python 3.8 for example)
# - install numpy: pip install numpy
# - inside boost dir:
#    ./bootstrap.sh --with-libraries=python
#    Check ./project-config.jam and edit it if needed as it should include the path to the right python install
#    ./b2
#       If you have trouble using the right python version or just pointing to one python lib
#       ("fatal error: pyconfig.h: No such file or directory")
#       you may instead try : CPLUS_INCLUDE_PATH=/somepath/anaconda3/include/pythonX.Y ./b2
#    ./b2 install
find_package(Boost 1.75.0 COMPONENTS python38 numpy38 REQUIRED)
# Python library, will produce a icm.so file
python_add_module(icm src/python-bridge.cpp
        src/icm.cpp include/mtt-icm/icm.hpp
        src/monte-carlo-icm.cpp include/mtt-icm/monte-carlo-icm.hpp
        src/math.cpp src/math.hpp)

target_include_directories(icm PRIVATE ${PYTHON_INCLUDE_DIRS})
target_include_directories(icm PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(icm PRIVATE "${mtt_icm_SOURCE_DIR}/include")

target_link_libraries(icm ${PYTHON_LIBRARY})
target_link_libraries(icm ${Boost_PYTHON_LIBRARY})
target_link_libraries(icm ${Boost_LIBRARIES})
target_link_libraries(icm wtl::sfmt)

#### Tests

add_subdirectory(test)

