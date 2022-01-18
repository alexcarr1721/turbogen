# CMake project file for turbogen

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(turbogen Fortran)

# Set the turbogen version
SET(VERSION 0.0.1)

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

# Include the module that sets fortran flags for DEBUG, TESTING, and RELEASE
# across many different compilers
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# Find MPI libraries
FIND_PACKAGE(MPI REQUIRED)

# Setup HDF5 libraries
SET(HDF5_USE_STATIC_LIBRARIES ON) # Link static libraries for performance
FIND_PACKAGE(HDF5 COMPONENTS Fortran)

# Setup FFTW
INCLUDE(${CMAKE_MODULE_PATH}/FindFFTW3.cmake)

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

