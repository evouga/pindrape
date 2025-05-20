# - Try to find the LIBSHELL library
# Once done this will define
#
#  LIBSHELL_FOUND - system has LIBSHELL 
#  LIBSHELL_INCLUDE_DIR - **the** LIBSHELL include directory
#  LIBSHELL_LIBRARIES - the LIBSHELL library directory
if(LIBSHELL_FOUND)
    return()
endif()

find_path(LIBSHELL_INCLUDE_DIR ElasticShell.h
    HINTS
        ENV LIBSHELL_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/libshell
        ${CMAKE_SOURCE_DIR}/../libshell
        ${CMAKE_SOURCE_DIR}/../../libshell
        /usr
        /usr/local        
    PATH_SUFFIXES include
)

find_library(LIBSHELL_LIB_DIR libshell
    HINTS
        ENV LIBSHELL_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/libshell
        ${CMAKE_SOURCE_DIR}/../libshell
        ${CMAKE_SOURCE_DIR}/../../libshell
        /usr
        /usr/local        
    PATH_SUFFIXES lib 
)

find_library(LIBSHELL_OPTIMIZATION_DIR optimization
    HINTS
        ENV LIBSHELL_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/libshell
        ${CMAKE_SOURCE_DIR}/../libshell
        ${CMAKE_SOURCE_DIR}/../../libshell
        /usr
        /usr/local        
    PATH_SUFFIXES lib 
)

set(LIBSHELL_LIBRARIES ${LIBSHELL_OPTIMIZATION_DIR} ${LIBSHELL_LIB_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libshell
    "\nlibshell not found"
    LIBSHELL_INCLUDE_DIR LIBSHELL_LIBRARIES)
mark_as_advanced(LIBSHELL_INCLUDE_DIR LIBSHELL_LIBRARIES)

