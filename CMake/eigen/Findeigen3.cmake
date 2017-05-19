# - Try to find the chull library.
# Input variables:
#     EIGEN_ROOT - chull installation root
#
# Output variables:
#     EIGEN_INCLUDE_DIR - include directories
#     EIGEN_FOUND - is 'YES' if library found, 'NO' otherwise

set(EIGEN_FOUND YES)

find_path(
    EIGEN_INCLUDE_DIR
    NAMES
    "eigen3/Eigen"
    HINTS
    "${EIGEN_ROOT}"
    ENV
    "EIGEN_ROOT"
    PATH_SUFFIXES
    "include"
)

if(NOT EIGEN_INCLUDE_DIR)
  set(EIGEN_FOUND NO)
else()
    message(STATUS "Found eigen: ${EIGEN_INCLUDE_DIR}")
endif()

# No need to find other libraries yet.



# find_path(EIGEN3_MODULE_DIR
#   NAMES eigen3
#   HINTS /opt/local/include ${EIGEN3_MODULE_SRC_DIR}/.. ../../
#   DOC "Location of the eigen3 root directory"
#   NO_DEFAULT_PATH)
#
# if( ${EIGEN3_MODULE_DIR} STREQUAL "EIGEN3_MODULE_DIR-NOTFOUND" )
#   message( FATAL_ERROR "Can't find eigen3! Specify with -DEIGEN3_MODULE_DIR" )
# endif()
#
# find_path(EIGEN3_MODULE_SRC_DIR
#   NAMES Eigen
#   HINTS ${EIGEN3_MODULE_DIR}
#   DOC "Location of eigen3's Eigen/ directory"
#   NO_DEFAULT_PATH)
#
# set(CMAKE_MODULE_PATH ${EIGEN3_MODULE_DIR} ${CMAKE_MODULE_PATH})
#
# # Configure eigen #
# include_directories(${EIGEN3_MODULE_DIR})
# # end #
