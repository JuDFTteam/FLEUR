cmake_minimum_required(VERSION 3.0)
project(FLEUR LANGUAGES Fortran)

include("../cmake/Generate_Schema.cmake")

#set(FLEURINPUT_COMPILEOPTS "-r8")

add_library(fleurinput	 STATIC
fleurinput_read_xml.f90    fleurinput_mpi_bc.f90
types_fleurinput_base.f90  types_input.f90      types_oneD.f90
types_coreSpecInput.f90    types_hybinp.f90     types_noco.f90
types_vacuum.f90           types_banddos.f90    types_wannier.f90
types_cell.f90             types_kpts.f90       types_fleurinput.f90
types_sliceplot.f90        types_sym.f90        types_atoms.F90
types_econfig.F90          types_field.F90      types_xcpot.F90
types_enparaXML.f90        types_forcetheo_data.f90
types_mpinp.f90 
types_xml.f90
calculator.f constants.f90 mpi_bc_tool.F90 ../io/xml/inputSchema.h ../io/xml/dropInputSchema.c ../io/xml/xmlInterfaceWrapper.c
)
#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/fleurinput")

if (DEFINED FLEURINPUT_COMPILEOPTS )
  target_compile_options(fleurinput PRIVATE ${FLEURINPUT_COMPILEOPTS})
endif()

target_link_libraries(fleurinput juDFT)
set_target_properties(fleurinput PROPERTIES Fortran_MODULE_DIRECTORY modules/fleurinput)
