cmake_minimum_required(VERSION 3.0)
project(FLEUR LANGUAGES Fortran)

add_library(fleurinput	 STATIC
 types_coreSpecInput.f90    types_hybrid.f90  types_noco.f90         types_vacuum.f90 fleurinput_read_xml.f90
types_banddos.f90        types_fleurinput_base.f90  types_input.f90   types_oneD.f90         types_wannier.f90
types_cell.f90           types_fleurinput.f90       types_kpts.f90    types_sliceplot.f90  types_sym.f90
types_atoms.F90  types_econfig.F90    types_field.F90    types_xcpot.F90 types_enparaXML.f90 types_forcetheo_data.f90 
types_xml.f90 calculator.f constants.f90 mpi_bc_tool.F90
)
#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/fleurinput")
target_link_libraries(fleurinput juDFT)
