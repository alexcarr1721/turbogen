# Find the FFTW3 library and include files

# Check the following environment variables, if they 
# are defined then assume they point to the top directory
# of FFTW
if( NOT FFTW_TOP_DIR AND DEFINED ENV{FFTWROOT} )
    set( FFTW_TOP_DIR $ENV{FFTWROOT} )
elseif( NOT FFTW_TOP_DIR AND DEFINED ENV{FFTW_ROOT} )
    set( FFTW_TOP_DIR $ENV{FFTW_ROOT} )
elseif( NOT FFTW_TOP_DIR AND DEFINED ENV{FFTWDIR} )
    set( FFTW_TOP_DIR $ENV{FFTWDIR} )
elseif( NOT FFTW_TOP_DIR AND DEFINED ENV{FFTW_DIR} )
    set( FFTW_TOP_DIR $ENV{FFTW_DIR} )
elseif( NOT FFTW_TOP_DIR AND DEFINED ENV{HPC_FFTW_DIR} )
    set( FFTW_TOP_DIR $ENV{HPC_FFTW_DIR} )
endif()

# Use PkgConfig as a last ditch effort
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT FFTW_TOP_DIR )
    pkg_check_modules( PKG_FFTW QUIET "fftw3f" )
endif()

# Set to find static libraries (usually these are the only libs
# available with FFTW)
set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )

if ( FFTW_TOP_DIR )
    find_library(FFTW3_LIB 
        NAMES "fftw3f" 
        PATHS ${FFTW_TOP_DIR}
        PATH_SUFFIXES "lib" "lib64"
        NO_DEFAULT_PATH
        )
    find_path(FFTW_INCLUDE_DIRS
        NAMES "fftw3.h"
        PATHS ${FFTW_TOP_DIR}
        PATH_SUFFIXES "include"
        NO_DEFAULT_PATH
        )
else()
    find_library(FFTW3_LIB
        NAMES "fftw3f"
        PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
        )
    find_path(FFTW_INCLUDE_DIRS
        NAMES "fftw3.h"
        PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
        )
endif()

# Add FFTW3_LIB to FFTW_LIBRARIES if found
if (FFTW3_LIB)
    set(FFTW3_LIB_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW3_LIB})
else()
    set(FFTW3_LIB_FOUND FALSE)
endif()

# 
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_INCLUDE_DIRS
    HANDLE_COMPONENTS
    )

mark_as_advanced(
    FFTW_INCLUDE_DIRS
    FFTW_LIBRARIES
    FFTW3_LIB
    )
