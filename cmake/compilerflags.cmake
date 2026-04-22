#cmake file to set compiler flags for some of the known compilers

if (NOT DEFINED FLEUR_STRICT_COMPILER_CHECK)
   set(FLEUR_STRICT_COMPILER_CHECK ON)
endif()

function(_fleur_report_compiler_issue reason details)
   set(_msg "Compiler check failed: ${reason}\n")
   string(APPEND _msg "Compiler ID: ${CMAKE_Fortran_COMPILER_ID}\n")
   string(APPEND _msg "Compiler version: ${CMAKE_Fortran_COMPILER_VERSION}\n")
   string(APPEND _msg "Details: ${details}")

   if (FLEUR_STRICT_COMPILER_CHECK)
      message(FATAL_ERROR "${_msg}")
   else()
      message(WARNING "${_msg}")
   endif()
endfunction()

set(_fleur_fc_perf "")
set(_fleur_fc_safe "-O2 -g")
set(_fleur_fc_debug "-O0 -g -DCPP_DEBUG")
set(_fleur_cc_perf "")
set(_fleur_cc_safe "-O2 -g")
set(_fleur_cc_debug "-O0 -g")

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -qmkl")
   set(_fleur_fc_perf "-O2 -g")
   set(_fleur_cc_perf "-O2 -g")

   if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "20.0.0.0")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -assume byterecl")
      set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_OLDINTEL")
      set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_OLDINTEL")
      _fleur_report_compiler_issue(
         "Intel compiler too old"
         "Version < 13.0.0.0 is unsupported for reliable FLEUR builds. Please upgrade your Intel compiler."
      )
   else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl  -assume byterecl -no-wrap-margin")
   endif()
   if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "19.0.0.0")
       set(_fleur_fc_debug "-C -traceback -O0 -g -DCPP_DEBUG -warn all")
   else()
       set(_fleur_fc_debug "-check assume,bounds,contiguous,format,output_conversion,stack -traceback -O0 -g -DCPP_DEBUG")
   endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC ")
   set(FLEUR_COMPILE_OPTIONS -Mlre -Mautoinline -Mvect=simd -Mcache_align -Mflushz)

   set(_fleur_fc_perf "-O3 -g")
   set(_fleur_cc_perf "-O3 -g")
   set(_fleur_fc_debug "-C -traceback -O0 -g -Mchkstk -gpu=debug -DCPP_DEBUG")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
   set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") #fix problem in cmake
   #CPU
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
   set(FLEUR_COMPILE_OPTIONS -mavx2 -Mlre -Mautoinline -Mpre -Mvect=simd -Mcache_align -Mflushz)

   set(_fleur_fc_perf "-O1 -g")
   set(_fleur_cc_perf "-O1 -g")
   set(_fleur_fc_debug "-C -traceback -O0 -g -Mchkstk -Mchkptr -Ktrap=fp -DCPP_DEBUG")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnosave -qarch=qp -qtune=qp -qfixed -qsuppress=1520-022 -qessl")

   set(_fleur_fc_perf "-O4 -qsuppress=1500-036")
   set(_fleur_cc_perf "-O3 -g")
   set(_fleur_fc_debug "-O0 -g -DCPP_DEBUG")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/bgsys/local/libxml2/include/libxml2")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_AIX")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_AIX")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "13.0")
      _fleur_report_compiler_issue(
         "GNU compiler too old"
         "Only modern versions of gfortran >=13 are supported for FLEUR. Please use a newer GNU toolchain."
      )
   endif()
   #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wno-missing-include-dirs -DCPP_IRAPPROX")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wno-missing-include-dirs -fno-sign-zero -fPIC")

   set(_fleur_fc_perf "-Ofast -g")
   set(_fleur_cc_perf "-Ofast -g")
   set(_fleur_fc_debug "-fdump-core -Wall -Wextra -Wno-array-temporaries -fbacktrace -fcheck=all -finit-real=nan -O0 -g -DCPP_DEBUG")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
 elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
   set(_fleur_fc_perf "-O3")
   set(_fleur_cc_perf "-O3")
   set(_fleur_fc_debug "-O0 -g -DCPP_DEBUG")
  elsif (CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -kind=byte")
   set(_fleur_fc_perf "-O2 -g")
   set(_fleur_cc_perf "-O2 -g")
   set(_fleur_fc_debug "-C=all -gline -O0 -DCPP_DEBUG")
  else()
   _fleur_report_compiler_issue(
      "Unsupported compiler"
      "No compiler flag mapping is implemented for this Fortran compiler ID."
   )
endif()

if (FLEUR_OPT_LEVEL STREQUAL "debug")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_fleur_fc_debug}")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_fleur_cc_debug}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_fleur_cc_debug}")
elseif (FLEUR_OPT_LEVEL STREQUAL "safe")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_fleur_fc_safe}")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_fleur_cc_safe}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_fleur_cc_safe}")
else()
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_fleur_fc_perf}")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_fleur_cc_perf}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_fleur_cc_perf}")
endif()

if (DEFINED FLEUR_FORTRAN_ARCH_FLAGS_SELECTED AND NOT "${FLEUR_FORTRAN_ARCH_FLAGS_SELECTED}" STREQUAL "")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLEUR_FORTRAN_ARCH_FLAGS_SELECTED}")
endif()
if (DEFINED FLEUR_C_ARCH_FLAGS_SELECTED AND NOT "${FLEUR_C_ARCH_FLAGS_SELECTED}" STREQUAL "")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FLEUR_C_ARCH_FLAGS_SELECTED}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLEUR_C_ARCH_FLAGS_SELECTED}")
endif()

