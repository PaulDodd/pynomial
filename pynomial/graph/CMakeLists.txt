set(_graph_sources   _graph.cxx
                    )

set(_graph_headers
    module.h
    Matrix.h
    Network.h
    )


# Need to define NO_IMPORT_ARRAY in every file but module.cc
set_source_files_properties(${_graph_sources} PROPERTIES COMPILE_DEFINITIONS NO_IMPORT_ARRAY)

# if (ENABLE_CUDA)
# CUDA_COMPILE(_CUDA_GENERATED_FILES ${_graph_cu_sources} OPTIONS ${CUDA_ADDITIONAL_OPTIONS} SHARED)
# endif (ENABLE_CUDA)

add_library (_graph SHARED ${_graph_sources})

# link the library to its dependencies
target_link_libraries(_graph _pynomial ${PYNOMIAL_COMMON_LIBS})

# if we are compiling with MPI support built in, set appropriate
# compiler/linker flags
if (ENABLE_MPI)
   if(MPI_COMPILE_FLAGS)
       set_target_properties(_graph PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
   endif(MPI_COMPILE_FLAGS)
   if(MPI_LINK_FLAGS)
       set_target_properties(_graph PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
   endif(MPI_LINK_FLAGS)
endif(ENABLE_MPI)

##################################
# change the name of the library to be a valid python module
# tweak the properties of the output to make a functional python module
set_target_properties(_graph PROPERTIES PREFIX "" OUTPUT_NAME "_graph")

# .dylib is not recognized as a python module by python on Mac OS X
if(APPLE)
    set_target_properties(_graph PROPERTIES SUFFIX ".so")
endif(APPLE)
fix_cudart_rpath(_graph)

# install the library
install(TARGETS _graph
        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/graph
        )

################ Python only modules
# copy python modules to the build directory to make it a working python package
MACRO(copy_file file)
    add_custom_command (
        OUTPUT ${file}
        DEPENDS ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}
        COMMENT    "Copy pynomial/graph/${file}"
    )
ENDMACRO(copy_file)

set(files  __init__.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/graph
       )

foreach(file ${files})
    copy_file(${file})
endforeach()

add_custom_target(copy_graph ALL DEPENDS ${files})

# copy python modules to the build directory to make it a working python package
MACRO(copy_header file)
    add_custom_command (
        OUTPUT ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/../include/pynomial/graph/${file}
        COMMENT    "Copy pynomial/include/pynomial/graph/${file}"
    )
ENDMACRO(copy_header)

# install headers in installation target
install(FILES ${_graph_headers}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/include/pynomial/graph
       )

if (COPY_HEADERS)
# copy headers into build directory
foreach(file ${_graph_headers})
    copy_header(${file})
endforeach()

add_custom_target(copy_graph_headers ALL DEPENDS ${_graph_headers})
endif()

if (BUILD_TESTING)
    add_subdirectory(test-py)
    add_subdirectory(test)
endif()

if (BUILD_VALIDATION)
    add_subdirectory(validation)
endif()
