########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(main_src ${SRCMAIN}/main2.f90)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${MAINEXE} ${main_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links main to the base library
TARGET_LINK_LIBRARIES(${MAINEXE} ${TURBOBASELIB})

# Uncomment if you need to link to BLAS and LAPACK
#TARGET_LINK_LIBRARIES(${FOOEXE} ${BLAS_LIBRARIES}
#                                ${LAPACK_LIBRARIES}
#                                ${CMAKE_THREAD_LIBS_INIT})

# Link HDF5 and MKL Libraries?
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${MAINEXE} ${HDF5_Fortran_LIBRARIES})

IF(USE_MPI)
    SET_TARGET_PROPERTIES(${MAINEXE} PROPERTIES
                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${MAINEXE} ${MPI_Fortran_LIBRARIES})
ENDIF(USE_MPI)

#####################################
# Tell how to install this executable
#####################################

INSTALL(TARGETS ${MAINEXE} DESTINATION bin)
