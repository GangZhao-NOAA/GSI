macro (setGeneric)
  message("Setting paths for Generic System")
  option(FIND_HDF5 "Try to Find HDF5 libraries" OFF)
  option(FIND_HDF5_HL "Try to Find HDF5 libraries" OFF)
  if(EXISTS /jetmon) 
    set(HOST_FLAG "" CACHE INTERNAL "Host Flag") ## default, no host_flag required
  else()
    set(HOST_FLAG "" CACHE INTERNAL "Host Flag") ## default, no host_flag required
  endif()
  set(MKL_FLAG "-mkl"  CACHE INTERNAL "MKL Flag")
  set(GSI_Intel_Platform_FLAGS "-DPOUND_FOR_STRINGIFY -O3 -fp-model source -assume byterecl -convert big_endian -g -traceback -D_REAL8_ ${OpenMP_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}" CACHE INTERNAL "GSI Fortran Flags")
  set(ENKF_Platform_FLAGS "-O3 ${HOST_FLAG} -warn all -implicitnone -traceback -fp-model strict -convert big_endian -DGFS -D_REAL8_  ${OpenMP_Fortran_FLAGS}" CACHE INTERNAL "ENKF Fortran Flags")

  message("setting values for corelibs")
  set(BUILD_BUFR "OFF" CACHE INTERNAL "Build the BUFR library" )
  set(BUILD_BACIO "OFF" CACHE INTERNAL "Build the BACIO library" )
  set(BUILD_SFCIO "OFF" CACHE INTERNAL "Build the SFCIO library" )
  set(BUILD_SIGIO "OFF" CACHE INTERNAL "Build the SIGIO library" )
  set(BUILD_NEMSIO "OFF" CACHE INTERNAL "Build the NEMSIO library" )
  set(BUILD_SP "OFF" CACHE INTERNAL "Build the SP library" )
  set(BUILD_CRTM "OFF" CACHE INTERNAL "Build the CRTM library" )
  set(BUILD_W3EMC "OFF" CACHE INTERNAL "Build the EMC library" )
  set(BUILD_W3NCO "OFF" CACHE INTERNAL "Build the EMC library" )
  set(BUILD_NCO "OFF" CACHE INTERNAL "Build the NCO library" )
endmacro()
