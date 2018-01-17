# configure_file (PYNOMIALVersion.h.in ${PYNOMIAL_BINARY_DIR}/pynomial/include/PYNOMIALVersion.h)

# install(FILES ${PYNOMIAL_BINARY_DIR}/pynomial/include/PYNOMIALVersion.h
#         DESTINATION ${PYTHON_MODULE_BASE_DIR}/include
#         )

set(_pynomial_sources extern/num_util.cc
                   )

set(_pynomial_headers
    )


# Need to define NO_IMPORT_ARRAY in every file but module.cc
set_source_files_properties(${_pynomial_sources} PROPERTIES COMPILE_DEFINITIONS NO_IMPORT_ARRAY)

add_library (_pynomial SHARED _pynomial.cxx ${_pynomial_sources})

# link the library to its dependencies
target_link_libraries(_pynomial ${PYNOMIAL_COMMON_LIBS})

set(PYNOMIAL_LIBRARIES _pynomial ${PYNOMIAL_COMMON_LIBS})

# if we are compiling with MPI support built in, set appropriate
# compiler/linker flags
if (ENABLE_MPI)
   if(MPI_COMPILE_FLAGS)
       set_target_properties(_pynomial PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
   endif(MPI_COMPILE_FLAGS)
   if(MPI_LINK_FLAGS)
       set_target_properties(_pynomial PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
   endif(MPI_LINK_FLAGS)
endif(ENABLE_MPI)

##################################
# change the name of the library to be a valid python module
# tweak the properties of the output to make a functional python module
set_target_properties(_pynomial PROPERTIES PREFIX "" OUTPUT_NAME "_pynomial")

# .dylib is not recognized as a python module by python on Mac OS X
if(APPLE)
    set_target_properties(_pynomial PROPERTIES SUFFIX ".so")
endif(APPLE)
fix_cudart_rpath(_pynomial)

# install the library
install(TARGETS _pynomial
        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}
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
        COMMENT    "Copy pynomial/${file}"
    )
ENDMACRO(copy_file)

set(files __init__.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}
       )

foreach(file ${files})
    copy_file(${file})
endforeach()

# copy python modules to the build directory to make it a working python package
MACRO(copy_header file)
    add_custom_command (
        OUTPUT ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/include/pynomial/${file}
        COMMENT    "Copy pynomial/include/pynomial/${file}"
    )
ENDMACRO(copy_header)

add_custom_target(copy_pynomial ALL DEPENDS ${files})

# install headers in installation target
install(FILES ${_pynomial_headers}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/include/pynomial
       )

if (COPY_HEADERS)
# copy headers into build directory
foreach(file ${_pynomial_headers})
    copy_header(${file})
endforeach()

add_custom_target(copy_pynomial_headers ALL DEPENDS ${_pynomial_headers})
endif()

file(GLOB _directory_contents RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

# explicitly remove packages which are already explicitly dealt with
# list(REMOVE_ITEM _directory_contents test)

foreach(entry ${_directory_contents})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${entry} OR IS_SYMLINK ${CMAKE_CURRENT_SOURCE_DIR}/${entry})
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${entry}/CMakeLists.txt)
            message("adding ${entry}")
            add_subdirectory(${entry})
        endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${entry}/CMakeLists.txt)
    endif(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${entry} OR IS_SYMLINK ${CMAKE_CURRENT_SOURCE_DIR}/${entry})
endforeach()
