if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   message("Intel Fortran detected")
   set(FLEUR_PRECISION_OPTION "-r8")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  message("NVHPC Fortran detected")
  set(FLEUR_PRECISION_OPTION "-Mr8;-Mr8intrinsics")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
   message("PGI Fortran detected")
   set(FLEUR_PRECISION_OPTION "-Mr8;-Mr8intrinsics")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
   message("IBM/BG Fortran detected")
   set(FLEUR_PRECISION_OPTION "-qrealsize=8")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   message("gfortran detected")
   set(FLEUR_PRECISION_OPTION "-fdefault-real-8;-fdefault-double-8")
 elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
	 message("Cray compiler detected")
	 set(FLEUR_PRECISION_OPTION -s real64)
else()
   message("Unknown compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
endif()

string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLEUR_PRECISION_OPTION}")
