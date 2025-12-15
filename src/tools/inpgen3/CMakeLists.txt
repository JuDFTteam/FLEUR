cmake_minimum_required(VERSION 3.5)
project(FLEUR LANGUAGES Fortran CXX)

#These options should be adjusted
option(INPGEN_USE_DEBUG "Set CPP_DEBUG when compiling code" OFF)
#In addition you might want to set
#set(JUDFT_LIBRARIES "-lsomething")
#set(JUFDT_COMPILEOPTS "-r8")
option(INPGEN_USE_kplib "Use kplib library" OFF)

set(FLEUR_SRC "${CMAKE_SOURCE_DIR}/src/fleur")

if (INPGEN_USE_kplib)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   set(externallibs
	kpoints/kplib_interface.F90
	${FLEUR_SRC}/external/kplib/fleur_interface.cpp
)
else()
   set(externallibs kpoints/kplib_interface_dummy.F90)
endif()


set(INPGEN3_SOURCES
Profiles/dropDefaults.c
Profiles/dropDefault2EConfig.c
Profiles/dropOxidesValidationEConfig.c
Profiles/dropProfiles.c
${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h

inpgen_help.f90
inpgen_version.f90
bravais_symm.f90
make_crystal.f90
make_spacegroup.f90
check_mt_radii.f90
make_atom_groups.f90
make_magnetism.f90
make_defaults.f90
process_lattice_namelist.f90
atompar.F90
film_sym.f90
make_atomic_defaults.f90
read_inpgen_input.f90
super_check.f90
make_kpoints.f90
Profiles/types_Profile.f90
mod_symdata.f
inpgen_make_inp.F90
inpgen_make_kpts.F90

${CMAKE_BINARY_DIR}/include/buildinfo.h

${FLEUR_SRC}/main/fleur_dropxmlschema.f90
${FLEUR_SRC}/init/make_sym.f90
${FLEUR_SRC}/init/angles.f90
${FLEUR_SRC}/init/compile_descr.F90
${FLEUR_SRC}/init/ss_sym.f
${FLEUR_SRC}/global/sort.f90
${FLEUR_SRC}/math/inv3.f90
w_inpXML.f90

kpoints/divi.f
#kpoints/brzone2.f90
kpoints/kpttet.f
kpoints/kptmop.f
kpoints/ordstar.f
kpoints/fulstar.f
#kpoints/tetcon.f90
kpoints/kvecon.f
kpoints/kprep.f
#kpoints/triang.f
kpoints/make_tetra.f

${FLEUR_SRC}/math/util.F
${FLEUR_SRC}/math/pol_angle.f90
${FLEUR_SRC}/init/spgrot.f90
${FLEUR_SRC}/init/mapatom.F90


${FLEUR_SRC}/types/types_nococonv.F90
${FLEUR_SRC}/types/types_mat.F90
${FLEUR_SRC}/types/types_xcpot_inbuild_nofunction.F90
${FLEUR_SRC}/types/types_xcpot_data.F90
${FLEUR_SRC}/types/types_xcpot_libxc.F90
${FLEUR_SRC}/types/types_mpi.F90
${FLEUR_SRC}/types/types_forcetheo.F90
${FLEUR_SRC}/types/types_stars.f90
${FLEUR_SRC}/types/types_lapw.F90
${FLEUR_SRC}/types/types_sphhar.f90
${FLEUR_SRC}/types/types_tlmplm.F90
${FLEUR_SRC}/types/types_misc.F90
${FLEUR_SRC}/types/types_mpimat.F90
${FLEUR_SRC}/types/types_parallelLoop.F90
${FLEUR_SRC}/types/types_potden.F90
${FLEUR_SRC}/types/types_enpara.F90
${FLEUR_SRC}/types/types_setup.F90
${FLEUR_SRC}/types/types_hub1data.f90

${FLEUR_SRC}/hybrid/divide_most_evenly.f90
kpoints/kptgen_hybrid.f90
${FLEUR_SRC}/io/nocoInputCheck.F90
${FLEUR_SRC}/eigen/orthoglo.F90
${FLEUR_SRC}/math/ylm4.F90
${FLEUR_SRC}/cdn/slab_dim.f90
${FLEUR_SRC}/cdn/slabgeom.f90
${FLEUR_SRC}/propcalc/dos/nstm3.f90
${FLEUR_SRC}/eigen/int_21.f90
${FLEUR_SRC}/eigen/int_21lo.f90

${FLEUR_SRC}/xc-pot/gaunt.f90
#${FLEUR_SRC}/math/grule.f90
${FLEUR_SRC}/init/ifft235.f90
kpoints/gkptwgt.f90
${FLEUR_SRC}/init/local_sym.f
${FLEUR_SRC}/init/lhcal.f
${FLEUR_SRC}/init/ptsym.f90
${FLEUR_SRC}/math/gaussp.f
${FLEUR_SRC}/init/gtest.f
${FLEUR_SRC}/io/xmlOutput.F90
${FLEUR_SRC}/global/utility.F90
${FLEUR_SRC}/global/radsra.f
${FLEUR_SRC}/global/differ.f
${FLEUR_SRC}/global/soc_sym.f
${FLEUR_SRC}/math/inwint.f
${FLEUR_SRC}/math/outint.f
${FLEUR_SRC}/math/intgr.F90
#${FLEUR_SRC}/math/d_wigner.F90
${FLEUR_SRC}/math/trapz.f90
${FLEUR_SRC}/global/find_enpara.f90
${FLEUR_SRC}/io/xsf_io.f90
${FLEUR_SRC}/init/boxdim.f90
${FLEUR_SRC}/mpi/mpi_reduce_tool.F90
${FLEUR_SRC}/../libraries/fleurinput/mpi_bc_tool.F90
)


add_executable(inpgen3 ${externallibs}
  inpgen.f90
  ${INPGEN3_SOURCES}
) 

add_library(fleurinpgen SHARED ${externallibs}
  inpgen_lib.F90
  ${INPGEN3_SOURCES}
) 

set_source_files_properties(${CMAKE_BINARY_DIR}/include/buildinfo.h PROPERTIES GENERATED TRUE)

target_compile_options(inpgen3 PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
target_compile_options(fleurinpgen PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
  
include_directories("${CMAKE_CURRENT_BINARY_DIR}")
#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/inpgen3")
set_target_properties(inpgen3 PROPERTIES Fortran_MODULE_DIRECTORY modules/inpgen3)
set_target_properties(inpgen3 PROPERTIES RUNTIME_OUTPUT_NAME inpgen3)
set_target_properties(inpgen3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_compile_definitions(inpgen3 PUBLIC ${FLEUR_DEFINITIONS} "CPP_INPGEN")
target_link_libraries(inpgen3 fleurinput juDFT)
target_link_libraries(inpgen3 ${FLEUR_LIBRARIES})
target_link_libraries(inpgen3 ${FLEUR_LINK_LIBRARIES})

#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/fleurinpgen")
set_target_properties(fleurinpgen PROPERTIES Fortran_MODULE_DIRECTORY modules/fleurinpgen)
set_target_properties(fleurinpgen PROPERTIES RUNTIME_OUTPUT_NAME fleurinpgen)
set_target_properties(fleurinpgen PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
target_compile_definitions(fleurinpgen PUBLIC ${FLEUR_DEFINITIONS} "CPP_INPGEN")
target_link_libraries(fleurinpgen fleurinput juDFT)
target_link_libraries(fleurinpgen ${FLEUR_LIBRARIES})
target_link_libraries(fleurinpgen ${FLEUR_LINK_LIBRARIES})

if (INPGEN_USE_kplib)
  target_link_libraries(inpgen3 kplib_static)
endif()
set_target_properties(inpgen3 PROPERTIES LINKER_LANGUAGE Fortran)

if (INPGEN_USE_kplib)
  target_link_libraries(fleurinpgen kplib_static)
endif()
set_target_properties(fleurinpgen PROPERTIES LINKER_LANGUAGE Fortran)

file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/default_econfig.h)
file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/default2_econfig.h)
file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/oxides_validation_econfig.h)
file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/profileConfig.h)
find_program(XXD_PROG xxd)

if (XXD_PROG)
  
    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    COMMAND ${XXD_PROG} -i default.econfig ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/default.econfig
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/
    COMMENT "Putting current default.econfig into default_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    COMMAND ${XXD_PROG} -i default2.econfig ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/default2.econfig
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/
    COMMENT "Putting current default2.econfig into default2_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    COMMAND ${XXD_PROG} -i oxides_validation.econfig ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/oxides_validation.econfig
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/
    COMMENT "Putting current oxides_validation.econfig into oxides_validation_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    COMMAND ${XXD_PROG} -i profile.config ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/profile.config
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/ 
    COMMENT "Putting current profile.config into profileConfig.h")
else()
    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/default_econfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of default_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/default2_econfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of default2_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/oxides_validation_econfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of oxides_validation_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen3/Profiles/profileConfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of profileConfig.h")

endif()
