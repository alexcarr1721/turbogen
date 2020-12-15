# CMake project file for SBABL

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(TURBOGEN Fortran)

# Set the turbogen version
SET(VERSION 0.1)

# Add local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set default installation directory
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
IF(WIN32)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}\\install" CACHE PATH "Set default installation path." FORCE)
ELSE()
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Set default installation path." FORCE)
ENDIF(WIN32)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Uncomment if it is required that Fortran 90 is supported
#IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
#ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(MKL_STATIC "Use static linking for MKL libraries" ON)
OPTION(MKL_NEEDEXTRA "Import extra libraries, scalapack and cdft" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# Setup HDF5 libraries
FIND_PACKAGE(HDF5 COMPONENTS Fortran)

# Setup MKL_DFTI
INCLUDE(${CMAKE_MODULE_PATH}/FindMKL.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(MAINEXE turbogen)

# Define the base library name
SET(TURBOBASELIB base)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(INC ${CMAKE_SOURCE_DIR}/include)
SET(SRCMAIN ${SRC}/main)
SET(SRCBASE ${SRC}/base)

# Have the .mod files placed in the include folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${INC})

# The source for the libraries and place them in the lib folder
ADD_SUBDIRECTORY(${SRCBASE} ${LIB})

# The source for the sbABL binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRCMAIN} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)

# Add an uninstall target to the Makefile
IF(NOT TARGET uninstall)
  CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
ENDIF(NOT TARGET uninstall)


# Install to install/ or a directory pointed to by CMAKE_INSTALL_DIR
INSTALL(TARGETS turbogen DESTINATION bin)
INSTALL(FILES "${CMAKE_SOURCE_DIR}/lib/libbase.a" DESTINATION lib)
INSTALL(DIRECTORY  "${CMAKE_SOURCE_DIR}/include" DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.mod")
INSTALL(DIRECTORY  "${CMAKE_SOURCE_DIR}/include" DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.smod")

