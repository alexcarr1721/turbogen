#########################################################################################
#											#
#	Cmake script to locate Intel Math Kernel Library. Adapted from:			#
#	https://repository.prace-ri.eu/git/CodeVault/hpc-kernels/structured_grids	#
#											#
#########################################################################################


# Stage 1: Find the root directory

set(MKLROOT_PATH $ENV{MKLROOT})

if (NOT MKLROOT_PATH)
	# try to find at /opt/intel/mkl	
	if (EXISTS "/opt/intel/mkl")
		set(MKLROOT_PATH "/opt/intel/mkl")
	# try to find at K-cluster location
	elseif (EXISTS "/usr/local/pkgs-modules/intel_2019.1.144/mkl")
		set(MKLROOT_PATH "/usr/local/pkgs-modules/intel_2019.1.144/mkl")
	else ()
		MESSAGE(FATAL_ERROR "$MKLROOT is not set and cannot be found at /opt/intel/mkl")
	endif (EXISTS "/opt/intel/mkl")
endif (NOT MKLROOT_PATH)


# Stage 2: Find include path and libraries
	
if (MKLROOT_PATH)
	# root-path found
	
	set(EXPECT_MKL_INCPATH "${MKLROOT_PATH}/include")
	
	if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
	    set(EXPECT_MKL_LIBPATH "${MKLROOT_PATH}/lib")
	endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
	
	set(EXPECT_ICC_LIBPATH "$ENV{ICC_LIBPATH}")
	
	if (CMAKE_SYSTEM_NAME MATCHES "Linux")
	    if (CMAKE_SIZEOF_VOID_P MATCHES 8)
	        set(EXPECT_MKL_LIBPATH "${MKLROOT_PATH}/lib/intel64")
	    else (CMAKE_SIZEOF_VOID_P MATCHES 8)
	        set(EXPECT_MKL_LIBPATH "${MKLROOT_PATH}/lib/ia32")
	    endif (CMAKE_SIZEOF_VOID_P MATCHES 8)
	endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
	
	# set include
	
	if (IS_DIRECTORY ${EXPECT_MKL_INCPATH})
	    set(MKL_INCLUDE_DIR ${EXPECT_MKL_INCPATH})
	endif (IS_DIRECTORY ${EXPECT_MKL_INCPATH})
	
	if (IS_DIRECTORY ${EXPECT_MKL_LIBPATH})
		set(MKL_LIBRARY_DIR ${EXPECT_MKL_LIBPATH})
	endif (IS_DIRECTORY ${EXPECT_MKL_LIBPATH})
	
	# find specific library files
		
	find_library(LIB_MKL_RT NAMES mkl_rt HINTS ${MKL_LIBRARY_DIR})
	find_library(LIB_PTHREAD NAMES pthread)	
	
endif (MKLROOT_PATH)

set(MKL_LIBRARY 
	${LIB_MKL_RT} 
	${LIB_PTHREAD})
	
# deal with QUIET and REQUIRED argument

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MKL DEFAULT_MSG 
    MKL_LIBRARY_DIR
    LIB_MKL_RT
    LIB_PTHREAD
    MKL_INCLUDE_DIR)
    
mark_as_advanced(LIB_MKL_RT LIB_PTHREAD MKL_INCLUDE_DIR)


