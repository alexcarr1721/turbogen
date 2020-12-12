#########################################################
#							#
#	Find the MPI Libraries when USE_MPI=ON, 	#
#	otherwise make sure the libraries are		#
#	not used when compiling code. As of now		#
# 	compilation without MPI fails.			#
#							#
#########################################################

IF (USE_MPI)
    # Find MPI
    IF (NOT MPI_Fortran_FOUND)
        FIND_PACKAGE (MPI REQUIRED)
    ENDIF (NOT MPI_Fortran_FOUND)
ELSE ()
    # Turn off MPI (Just in case)
    UNSET (MPI_FOUND CACHE)
    UNSET (MPI_COMPILER CACHE)
    UNSET (MPI_LIBRARY CACHE)
ENDIF (USE_MPI)
