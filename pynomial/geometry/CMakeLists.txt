set(module_name "geometry")
set(_module_sources
                    )
set(_module_headers
    module.h
    sphere.h
    )

set(files   __init__.py
            sphere.py
    )


# Need to define NO_IMPORT_ARRAY in every file but module.cc
set_source_files_properties(${_module_sources} PROPERTIES COMPILE_DEFINITIONS NO_IMPORT_ARRAY)

add_library (_${module_name} SHARED _${module_name}.cxx)

# link the library to its dependencies
target_link_libraries(_${module_name} _pynomial ${PYNOMIAL_COMMON_LIBS})

# if we are compiling with MPI support built in, set appropriate
# compiler/linker flags
if (ENABLE_MPI)
   if(MPI_COMPILE_FLAGS)
       set_target_properties(_${module_name} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
   endif(MPI_COMPILE_FLAGS)
   if(MPI_LINK_FLAGS)
       set_target_properties(_${module_name} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
   endif(MPI_LINK_FLAGS)
endif(ENABLE_MPI)

##################################
# change the name of the library to be a valid python module
# tweak the properties of the output to make a functional python module
set_target_properties(_${module_name} PROPERTIES PREFIX "" OUTPUT_NAME "_${module_name}")

# .dylib is not recognized as a python module by python on Mac OS X
if(APPLE)
    set_target_properties(_${module_name} PROPERTIES SUFFIX ".so")
endif(APPLE)
fix_cudart_rpath(_${module_name})

# install the library
install(TARGETS _${module_name}
        LIBRARY DESTINATION ${PYTHON_MODULE_BASE_DIR}/${module_name}
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
        COMMENT    "Copy pynomial/${module_name}/${file}"
    )
ENDMACRO(copy_file)

install(FILES ${files}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/${module_name}
       )

foreach(file ${files})
    copy_file(${file})
endforeach()

add_custom_target(copy_${module_name} ALL DEPENDS ${files})

# copy python modules to the build directory to make it a working python package
MACRO(copy_header file)
    add_custom_command (
        OUTPUT ${file}
        POST_BUILD
        COMMAND    ${CMAKE_COMMAND}
        ARGS       -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/../include/pynomial/${module_name}/${file}
        COMMENT    "Copy pynomial/include/pynomial/${module_name}/${file}"
    )
ENDMACRO(copy_header)

# install headers in installation target
install(FILES ${_module_headers}
        DESTINATION ${PYTHON_MODULE_BASE_DIR}/include/pynomial/${module_name}
       )

if (COPY_HEADERS)
# copy headers into build directory
foreach(file ${_module_headers})
    copy_header(${file})
endforeach()

add_custom_target(copy_${module_name}_headers ALL DEPENDS ${_module_headers})
endif()

# if (BUILD_TESTING)
#     add_subdirectory(test)
# endif()
