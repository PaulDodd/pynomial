set(_polyhedron_sources
                    )

set(_polyhedron_headers
    module.h
    Polyhedra.h
    Net.h
    )


# Need to define NO_IMPORT_ARRAY in every file but module.cc
set_source_files_properties(${_polyhedron_sources} PROPERTIES COMPILE_DEFINITIONS NO_IMPORT_ARRAY)

# if (ENABLE_CUDA)
# CUDA_COMPILE(_CUDA_GENERATED_FILES ${_polyhedron_cu_sources} OPTIONS ${CUDA_ADDITIONAL_OPTIONS} SHARED)
# endif (ENABLE_CUDA)

add_library (_polyhedron SHARED _polyhedron.cxx)

# link the library to its dependencies
target_link_libraries(_polyhedron _graph _pynomial ${PYNOMIAL_COMMON_LIBS})

# if we are compiling with MPI support built in, set appropriate
# compiler/linker flags
if (ENABLE_MPI)
   if(MPI_COMPILE_FLAGS)
       set_target_properties(_polyhedron PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
   endif(MPI_COMPILE_FLAGS)
   if(MPI_LINK_FLAGS)
       set_target_properties(_polyhedron PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
   endif(MPI_LINK_FLAGS)
endif(ENABLE_MPI)

##################################
# change the name of the library to be a valid python module
# tweak the properties of the output to make a functional python module
set_target_properties(_polyhedron PROPERTIES PREFIX "" OUTPUT_NAME "_polyhedron")

# .dylib is not recognized as a python module by python on Mac OS X
if(APPLE)
    set_target_properties(_polyhedron PROPERTIES SUFFIX ".so")
endif(APPLE)
fix_cudart_rpath(_polyhedron)

# install the library
install(TARGETS _polyhedron
        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/polyhedron
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
        COMMENT    "Copy pynomial/polyhedron/${file}"
    )
ENDMACRO(copy_file)

set(files   __init__.py
            polyhedron.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/polyhedron
       )

foreach(file ${files})
    copy_file(${file})
endforeach()

add_custom_target(copy_polyhedron ALL DEPENDS ${files})

# copy python modules to the build directory to make it a working python package
MACRO(copy_header file)
    add_custom_command (
        OUTPUT ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/../include/pynomial/polyhedron/${file}
        COMMENT    "Copy pynomial/include/pynomial/polyhedron/${file}"
    )
ENDMACRO(copy_header)

# install headers in installation target
install(FILES ${_polyhedron_headers}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/include/pynomial/polyhedron
       )

if (COPY_HEADERS)
# copy headers into build directory
foreach(file ${_polyhedron_headers})
    copy_header(${file})
endforeach()

add_custom_target(copy_polyhedron_headers ALL DEPENDS ${_polyhedron_headers})
endif()

if (BUILD_TESTING)
    add_subdirectory(test-py)
    add_subdirectory(test)
endif()

if (BUILD_VALIDATION)
    add_subdirectory(validation)
endif()
