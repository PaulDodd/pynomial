cmake_minimum_required (VERSION 2.8)
if (CMAKE_MAJOR_VERSION VERSION_GREATER 2)
    if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
        cmake_policy(SET CMP0026 OLD)
        cmake_policy(SET CMP0042 OLD)
    endif(COMMAND cmake_policy)
endif()
project (pynomial)

#set some options
# option(PYNOMIAL_BUILD_SHARED_LIBS "Build shared libraries." ON)
# set(PYNOMIAL_BUILD_SHARED_LIBS ON) # right now this is the only option.
# option(PYNOMIAL_INSTALL_PYTHON_USER "Install the python module into the user site-packages directory" ON)

# The version number.
set (PYNOMIAL_VERSION_MAJOR 0)
set (PYNOMIAL_VERSION_MINOR 0)

# find dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

#find python
include(PythonSetup)
# add include directories
include_directories(${PYTHON_INCLUDE_DIRS})
list(APPEND LINK_LIBRARIES ${PYTHON_LIBRARIES})

# find eigen3
find_package(eigen3)
find_package(jansson)
if(JANSSON_FOUND)
    include_directories(${JANSSON_INCLUDE_DIRS})
    list(APPEND LINK_LIBRARIES ${JANSSON_LIBRARIES})
endif(JANSSON_FOUND)

# setup the install directories
set(DATA_INSTALL_DIR "share/${PROJECT_NAME}")
set(LIB_INSTALL_DIR "lib${LIB_SUFFIX}/${PROJECT_NAME}")
set(LIB_BASE_INSTALL_DIR "lib${LIB_SUFFIX}")
set(BIN_INSTALL_DIR "bin")
set(INC_INSTALL_DIR "include/${PROJECT_NAME}")
set(PYNOMIAL_PYTHON_MODULE_DIR ${PYTHON_SITEDIR}/${PROJECT_NAME})

set(ICP_SRC_DIR "extern/libicp/src")

# The project directories.
include_directories(${PROJECT_INCLUDE_DIR})

#end
