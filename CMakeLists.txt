cmake_minimum_required (VERSION 2.8)
if (CMAKE_MAJOR_VERSION VERSION_GREATER 2)
    if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
        cmake_policy(SET CMP0026 OLD)
        cmake_policy(SET CMP0042 OLD)
    endif(COMMAND cmake_policy)
endif()
project (procrustes)

#set some options
option(PROCRUSTES_BUILD_SHARED_LIBS "Build shared libraries." ON)
set(PROCRUSTES_BUILD_SHARED_LIBS ON) # right now this is the only option.
option(PROCRUSTES_INSTALL_PYTHON_USER "Install the python module into the user site-packages directory" ON)

# The version number.
set (PROCRUSTES_VERSION_MAJOR 0)
set (PROCRUSTES_VERSION_MINOR 0)

# find dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

#find python
include(PythonSetup)
# find boost
set(BOOST_PYTHON_COMPONENT "python")
set(REQUIRED_BOOST_COMPONENTS thread filesystem ${BOOST_PYTHON_COMPONENT}) # signals program_options unit_test_framework iostreams serialization
find_package(Boost 1.50.0 COMPONENTS REQUIRED ${REQUIRED_BOOST_COMPONENTS})
# add include directories
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
list(APPEND LINK_LIBRARIES ${Boost_LIBRARIES})
include_directories(${PYTHON_INCLUDE_DIRS})
list(APPEND LINK_LIBRARIES ${PYTHON_LIBRARIES})

# find eigen3
find_package(eigen3)
# find the folding project
find_package(sharedfolding)
# find the utilities
find_package(utilities)
# find the utilities
message("Finding jansson pre...")
find_package(jansson)
if(JANSSON_FOUND)
    include_directories(${JANSSON_INCLUDE_DIRS})
    list(APPEND LINK_LIBRARIES ${JANSSON_LIBRARIES})
endif(JANSSON_FOUND)
# find the utilities
find_package(json)


# setup the install directories
set(DATA_INSTALL_DIR "share/${PROJECT_NAME}")
set(LIB_INSTALL_DIR "lib${LIB_SUFFIX}/${PROJECT_NAME}")
set(LIB_BASE_INSTALL_DIR "lib${LIB_SUFFIX}")
set(BIN_INSTALL_DIR "bin")
set(INC_INSTALL_DIR "include/${PROJECT_NAME}")
set(PROCRUSTES_PYTHON_MODULE_DIR ${PYTHON_SITEDIR}/${PROJECT_NAME})

# message("Install directories: ${DATA_INSTALL_DIR} ${LIB_INSTALL_DIR} ${LIB_BASE_INSTALL_DIR} ${LIB_BASE_INSTALL_DIR} ${BIN_INSTALL_DIR} ${INC_INSTALL_DIR}")
# message("Python Install directories: ${PROCRUSTES_PYTHON_MODULE_DIR}")

set(ICP_SRC_DIR "extern/libicp/src")
add_subdirectory(cpp)
add_subdirectory(py)

# message("module path: ${CMAKE_MODULE_PATH}")

# The project directories.
include_directories(${PROJECT_INCLUDE_DIR})

# configure_file (
#   "${PROJECT_SOURCE_DIR}/FoldingConfig.h.in"
#   "${PROJECT_SOURCE_DIR}/FoldingConfig.h"
# )

# set(BUILD_TESTS OFF)
# if(BUILD_TESTS)
#     add_subdirectory(${PROJECT_TEST_SOURCE_DIR})
# endif(BUILD_TESTS)



















#end
