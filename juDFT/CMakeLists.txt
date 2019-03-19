


if (${FLEUR_USE_HDF5})
set(judft_HDF 
juDFT/hdf_accessprp.F90
juDFT/hdf_tools.F90
juDFT/hdf_tools_attr.F90
juDFT/hdf_tools_exists.F90
juDFT/hdf_tools_misc.F90
juDFT/hdf_tools_rw.F90
juDFT/hdf_tools_rw_var.F90
juDFT/hdf_tools_stride.F90)
else()
set(judft_HDF "")
endif()

add_library(juDFT STATIC
${judft_HDF} 
juDFT/usage_data.F90
juDFT/info.F90
juDFT/init.F90
juDFT/juDFT.F90
juDFT/stop.F90
juDFT/string.f90
juDFT/time.F90
juDFT/args.F90
juDFT/sysinfo.F90
juDFT/check_arguments.F90
juDFT/xmlOutput.F90
)
target_compile_definitions(juDFT PUBLIC ${FLEUR_DEFINITIONS})
target_link_libraries(juDFT ${FLEUR_LIBRARIES})
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/juDFT")
set_target_properties(juDFT PROPERTIES Fortran_MODULE_DIRECTORY modules/juDFT)
