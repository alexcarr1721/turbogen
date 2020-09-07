SHELL = /bin/sh

FC = ifort
CFLAGS  = -O2 -I${MKLROOT}/include -I$(LOCAL_HDF5_INC)
LFLAGS  =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a \
-Wl,--end-group -lpthread -L$(LOCAL_HDF5_LIB) \
$(LOCAL_HDF5_LIB)/libhdf5hl_fortran.a \
$(LOCAL_HDF5_LIB)/libhdf5_hl.a \
$(LOCAL_HDF5_LIB)/libhdf5_fortran.a \
$(LOCAL_HDF5_LIB)libhdf5.a -L$(LOCAL_ZLIB_LIB) \
-lz -ldl -lm -Wl,-rpath -Wl,$(LOCAL_HDF5_LIB)

default: rgm

rgm:  rgm.o functions.o mkl_dfti.o special_functions.o
	$(FC) $(CFLAGS) -o $@ $^ $(LFLAGS)

rgm.o: rgm.f90 functions.mod
	$(FC) $(CFLAGS) -c rgm.f90

functions.mod: functions.f90 functions.o
	@true

functions.o:  functions.f90 mkl_dfti.mod mkl_dft_type.mod special_functions.mod
	$(FC) $(CFLAGS) -c functions.f90

special_functions.mod: special_functions.f90 special_functions.o
	@true

special_functions.o: special_functions.f90
	$(FC) $(CFLAGS) -c special_functions.f90

mkl_dfti.mod: $(MKLROOT)/include/mkl_dfti.f90 mkl_dfti.o
	@true

mkl_dft_type.mod: $(MKLROOT)/include/mkl_dfti.f90 mkl_dfti.o
	@true

mkl_dfti.o:  $(MKLROOT)/include/mkl_dfti.f90
	$(FC) $(CFLAGS) -c $(MKLROOT)/include/mkl_dfti.f90

clean:
	$(RM) rgm *.o *~
