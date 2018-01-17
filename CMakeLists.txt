CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)

if (CMAKE_MAJOR_VERSION VERSION_GREATER 2)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
    cmake_policy(SET CMP0026 OLD)
    cmake_policy(SET CMP0042 OLD)
    set(CMAKE_MACOSX_RPATH ON)
endif(COMMAND cmake_policy)
endif()

project (PYNOMIAL)
option(INSTALL_SITE OFF "install pynomial as a site package")
# bring in custom modules
add_subdirectory (CMake)

#################################
## setup python library and executable
# setup MPI support
# include (MPISetup)
# find the python libraries to link to
include(PNPythonSetup)

if (INSTALL_SITE)
    SET(CMAKE_INSTALL_PREFIX ${PYTHON_SYSTEM_SITE} CACHE PATH "Python site installation directory" FORCE)
    message(STATUS "Setting installation site dir: ${CMAKE_INSTALL_PREFIX}")
ENDIF()

set(PYTHON_MODULE_BASE_DIR "pynomial")

# Find CUDA and set it up
# include (PNCUDASetup)
# Set default CFlags
include (PNCFlagsSetup)
# include some os specific options
# include (PNOSSpecificSetup)
# setup common libraries used by all targets in this project
include (PNCommonLibsSetup)
# setup numpy
include (PNNumpySetup)
# macros
include (PNMacros)

find_package(chull REQUIRED)
if(NOT CHULL_FOUND)
    message( FATAL_ERROR "Could not find chull." )
endif()
include_directories(${CHULL_INCLUDE_DIR} SYSTEM)
find_package(eigen3 REQUIRED)
if(NOT EIGEN_FOUND)
    message( FATAL_ERROR "Could not find eigen." )
endif()
include_directories(${EIGEN_INCLUDE_DIR} SYSTEM)

# Configure some source files, include directories, and create variables listing all source files
#################################
## Setup include directories and file lists for sub directories
include_directories(${PYNOMIAL_SOURCE_DIR}
                    ${PYNOMIAL_SOURCE_DIR}/pynomial/extern/cereal/include
                    ${PYNOMIAL_BINARY_DIR}/pynomial/include)

exec_program("date +%x" OUTPUT_VARIABLE COMPILE_DATE)

# preprocessor flags
# include (CMake_preprocessor_flags.cmake)

set (UPDATE_SUBMODULES 0 CACHE BOOL "Update git submodules for dependency projects")

if (UPDATE_SUBMODULES)
  message(STATUS "Updating git submodules")
  execute_process(COMMAND git submodule update --init
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
else (UPDATE_SUBMODULES)
  message(STATUS "Not updating git submodules")
endif (UPDATE_SUBMODULES)

set(cmake_files CMake/pynomial/PNCFlagsSetup.cmake
                CMake/pynomial/PNCommonLibsSetup.cmake
                CMake/pynomial/test.cc
                CMake/pynomial/PNMacros.cmake
                CMake/pynomial/PNNumpySetup.cmake
                CMake/pynomial/PNPythonSetup.cmake
   )

# copy cmake scripts we currently use into build directory
MACRO(copy_cmake file)
    add_custom_command (
        OUTPUT pynomial/${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/pynomial/${file}
        COMMENT    "Copy pynomial/${file}"
    )
ENDMACRO(copy_cmake)

# install headers in installation target
install(FILES ${cmake_files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}
       )

# copy headers into build directory
foreach(file ${cmake_files})
    copy_cmake(${file})
endforeach()

add_custom_target(copy_cmake ALL
    DEPENDS ${cmake_files})

################################
## Process subdirectories
add_subdirectory (pynomial)
