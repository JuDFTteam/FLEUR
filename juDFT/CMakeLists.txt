cmake_minimum_required(VERSION 3.0)
project(FLEUR LANGUAGES Fortran)

#These options should be adjusted
option(JUDFT_USE_MPI "Compile with MPI, will also work in serial" OFF)
option(JUDFT_USE_HDF5 "Compile HDF5 tools" OFF)
option(JUDFT_USE_HDF5MPI "Is the HDF5 version able to do parallel IO" OFF)
option(JUDFT_USE_COLLECT "Collect usage data" ON)
option(JUDFT_USE_DEBUG "Set CPP_DEBUG when compiling code" OFF)
#In addition you might want to set
#set(JUDFT_LIBRARIES "-lsomething")
#set(JUFDT_COMPILEOPTS "-r8")


if (${JUDFT_USE_HDF5})
set(judft_HDF 
hdf_accessprp.F90
hdf_tools.F90
hdf_tools_attr.F90
hdf_tools_exists.F90
hdf_tools_misc.F90
hdf_tools_rw.F90
hdf_tools_rw_var.F90
hdf_tools_stride.F90)
else()
set(judft_HDF "")
endif()

add_library(juDFT STATIC
${judft_HDF} 
usage_data.F90
info.F90
init.F90
juDFT.F90
stop.F90
string.f90
time.F90
args.F90
sysinfo.F90
check_arguments.F90
xmlOutput.F90
)
if (JUDFT_USE_MPI) 
  target_compile_definitions(juDFT PUBLIC CPP_MPI)
endif()
if (JUDFT_USE_HDF5)
  target_compile_definitions(juDFT PUBLIC CPP_HDF)
  target_link_libraries(juDFT "hdf5;hdf5_fortran")
endif()
if (JUDFT_USE_HDF5MPI)
  target_compile_definitions(juDFT PUBLIC CPP_HDFMPI)
endif()
if (JUDFT_USE_COLLECT)
  target_compile_definitions(juDFT PUBLIC CPP_ALLOW_USAGE_DATA)
endif()
if (JUDFT_USE_DEBUG)
  target_compile_definitions(juDFT PUBLIC CPP_DEBUG)
endif()

if (DEFINED JUDFT_LIBRARIES )
   target_link_libraries(juDFT ${JUDFT_LIBRARIES})
endif()

if (DEFINED JUDFT_COMPILEOPTS )
  target_compile_options(juDFT PUBLIC ${JUDFT_COMPILEOPTS})
endif()

#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/juDFT")
set_target_properties(juDFT PROPERTIES Fortran_MODULE_DIRECTORY modules/juDFT)
