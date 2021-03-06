include_directories(${PYNOMIAL_PYTHON_INCLUDE_DIR})

################################
## Define common libraries used by every target in PYNOMIAL

## to an ancient post online, adding -lutil fixed this in python 2.2
set(ADDITIONAL_LIBS "")
if (UNIX AND NOT APPLE)
    find_library(UTIL_LIB util /usr/lib)
    find_library(DL_LIB dl /usr/lib)
    set(ADDITIONAL_LIBS ${UTIL_LIB} ${DL_LIB})
    if (DL_LIB AND UTIL_LIB)
    mark_as_advanced(UTIL_LIB DL_LIB)
    endif (DL_LIB AND UTIL_LIB)
endif (UNIX AND NOT APPLE)

set(PYNOMIAL_COMMON_LIBS
        ${PYNOMIAL_PYTHON_LIBRARY}
        ${ADDITIONAL_LIBS}
        )

if (ENABLE_CUDA)
    list(APPEND PYNOMIAL_COMMON_LIBS ${CUDA_LIBRARIES} ${CUDA_cufft_LIBRARY} ${CUDA_curand_LIBRARY})

    if (ENABLE_NVTOOLS)
        list(APPEND PYNOMIAL_COMMON_LIBS ${CUDA_nvToolsExt_LIBRARY})
    endif()
endif (ENABLE_CUDA)

if (ENABLE_MPI)
    list(APPEND PYNOMIAL_COMMON_LIBS ${MPI_CXX_LIBRARIES})
endif (ENABLE_MPI)
