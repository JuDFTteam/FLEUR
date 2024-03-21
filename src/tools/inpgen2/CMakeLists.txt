cmake_minimum_required(VERSION 3.0)
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

add_executable(inpgen2 ${externallibs}
dropDefaults.c
dropDefault2EConfig.c
dropOxidesValidationEConfig.c
dropProfiles.c
default_econfig.h
default2_econfig.h
oxides_validation_econfig.h
profileConfig.h

inpgen_help.f90
inpgen_version.f90
inpgen.f90
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
read_old_inp.f90
types_Profile.f90
mod_symdata.f
old_inp/apws_dim.f90
old_inp/dimen7.F90
old_inp/dimens.F90
old_inp/first_glance.f
old_inp/fleur_init_old.F90
old_inp/inped.F90
old_inp/inpeig_dim.f90
old_inp/inpnoco.F90
old_inp/rw_inp.f90
old_inp/setlomap.F90
old_inp/setup.f90
old_inp/rw_noco.f90
old_inp/rw_symfile.f
old_inp/spg2set.f
old_inp/inpeig.f90
old_inp/strgn_dim.F
old_inp/prp_xcfft_box.f

${CMAKE_BINARY_DIR}/include/buildinfo.h

${FLEUR_SRC}/main/fleur_dropxmlschema.f90
${FLEUR_SRC}/init/make_sym.f90
${FLEUR_SRC}/init/angles.f90
${FLEUR_SRC}/init/compile_descr.F90
${FLEUR_SRC}/init/ss_sym.f
${FLEUR_SRC}/global/sort.f90
${FLEUR_SRC}/math/inv3.f90
${FLEUR_SRC}/io/w_inpXML.f90

kpoints/divi.f
#kpoints/bravais.f
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
${FLEUR_SRC}/cdn/int_21.f90
${FLEUR_SRC}/cdn/int_21lo.f90

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

set_source_files_properties(${CMAKE_BINARY_DIR}/include/buildinfo.h PROPERTIES GENERATED TRUE)

target_compile_options(inpgen2 PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
  
include_directories("${CMAKE_CURRENT_BINARY_DIR}")
#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/inpgen")
set_target_properties(inpgen2 PROPERTIES Fortran_MODULE_DIRECTORY modules/inpgen)
set_target_properties(inpgen2 PROPERTIES RUNTIME_OUTPUT_NAME inpgen)
set_target_properties(inpgen2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/..")
target_compile_definitions(inpgen2 PUBLIC ${FLEUR_DEFINITIONS} "CPP_INPGEN")
target_link_libraries(inpgen2 fleurinput juDFT)
target_link_libraries(inpgen2 ${FLEUR_LIBRARIES})
target_link_libraries(inpgen2 ${FLEUR_LINK_LIBRARIES})

if (INPGEN_USE_kplib)
  target_link_libraries(inpgen2 kplib_static)
endif()
set_target_properties(inpgen2 PROPERTIES LINKER_LANGUAGE Fortran)

file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/default_econfig.h)
file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/default2_econfig.h)
file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/oxides_validation_econfig.h)
file(REMOVE ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/profileConfig.h)
find_program(XXD_PROG xxd)

if (XXD_PROG)
    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    COMMAND ${XXD_PROG} -i default.econfig ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/default.econfig
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/
    COMMENT "Putting current default.econfig into default_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    COMMAND ${XXD_PROG} -i default2.econfig ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/default2.econfig
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/
    COMMENT "Putting current default2.econfig into default2_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    COMMAND ${XXD_PROG} -i oxides_validation.econfig ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/oxides_validation.econfig
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/
    COMMENT "Putting current oxides_validation.econfig into oxides_validation_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    COMMAND ${XXD_PROG} -i profile.config ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    DEPENDS ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/profile.config
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/
    COMMENT "Putting current profile.config into profileConfig.h")
else()
    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/default_econfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/default_econfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of default_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/default2_econfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/default2_econfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of oxides_validation_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/oxides_validation_econfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/oxides_validation_econfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of oxides_validation_econfig.h")

    ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    COMMAND cp ${CMAKE_SOURCE_DIR}/src/tools/inpgen2/profileConfig.h.backup ${CMAKE_CURRENT_BINARY_DIR}/profileConfig.h
    COMMENT "No xxd found using backup")
    message("No xxd command found! Using backup of profileConfig.h")

endif()
