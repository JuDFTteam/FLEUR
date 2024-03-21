cmake_minimum_required(VERSION 3.0)
project(FLEUR LANGUAGES Fortran)

include("${CMAKE_SOURCE_DIR}/cmake/Generate_Schema.cmake")

#set(FLEURINPUT_COMPILEOPTS "-r8")

set (fleurinp_src fleurinput_read_xml.f90    fleurinput_mpi_bc.f90
types_fleurinput_base.F90  types_input.f90      
types_coreSpecInput.f90    types_hybinp.f90     types_noco.f90
types_vacuum.f90           types_banddos.f90    types_wannier.f90
types_cell.f90             types_kpts.F90       types_fleurinput.f90
types_sliceplot.f90        types_sym.f90        types_atoms.F90
types_econfig.F90          types_field.F90      types_xcpot.F90
types_enparaXML.f90        types_forcetheo_data.f90
types_mpinp.f90            types_gfinp.F90      types_hub1inp.f90
types_xml.f90              types_juPhon.F90     types_brZone.f90
calculator.f90 constants.f90 mpi_bc_tool.F90)

set(fleurinp_src ${fleurinp_src}
${CMAKE_SOURCE_DIR}/src/fleur/io/xml/dropInputSchema.c 
${CMAKE_SOURCE_DIR}/src/fleur/io/xml/dropOutputSchema.c ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/inputSchema_old.h ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/outputSchema_old.h ${CMAKE_SOURCE_DIR}/src/fleur/io/xml/xmlInterfaceWrapper.c
${CMAKE_SOURCE_DIR}/src/fleur/math/d_wigner.F90 ${CMAKE_SOURCE_DIR}/src/fleur/math/inv3.f90 ${CMAKE_SOURCE_DIR}/src/fleur/math/grule.f90 ${CMAKE_SOURCE_DIR}/src/fleur/global/sort.f90 ${CMAKE_SOURCE_DIR}/src/fleur/global/atom_shells.f90 ${CMAKE_SOURCE_DIR}/src/fleur/kpoints/tetcon.f90 
${CMAKE_SOURCE_DIR}/src/fleur/kpoints/triang.f ${CMAKE_SOURCE_DIR}/src/fleur/kpoints/bravais.f ${CMAKE_SOURCE_DIR}/src/fleur/kpoints/brzone2.f90
${CMAKE_SOURCE_DIR}/src/fleur/hybrid/map_to_unit.f90
${CMAKE_BINARY_DIR}/include/inputSchema.h ${CMAKE_BINARY_DIR}/include/outputSchema.h 
)


add_library(fleurinput	 STATIC ${fleurinp_src})
#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/fleurinput")

if (DEFINED FLEURINPUT_COMPILEOPTS )
  target_compile_options(fleurinput PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEURINPUT_COMPILEOPTS}>)
endif()
target_compile_options(fleurinput PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)

target_compile_definitions(fleurinput PUBLIC ${FLEUR_DEFINITIONS})
target_link_libraries(fleurinput ${FLEUR_LINK_LIBRARIES})
target_link_libraries(fleurinput juDFT)
set_target_properties(fleurinput PROPERTIES Fortran_MODULE_DIRECTORY modules/fleurinput)

if ((CMAKE_Fortran_COMPILER_ID MATCHES "PGI") OR (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC"))
	set_source_files_properties(types_kpts.F90 TARGET_DIRECTORY fleurinput PROPERTIES COMPILE_OPTIONS "-nomp")
endif()    

