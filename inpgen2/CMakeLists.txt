cmake_minimum_required(VERSION 3.0)
project(FLEUR LANGUAGES Fortran)

#These options should be adjusted
option(INPGEN_USE_DEBUG "Set CPP_DEBUG when compiling code" OFF)
#In addition you might want to set
#set(JUDFT_LIBRARIES "-lsomething")
#set(JUFDT_COMPILEOPTS "-r8")


set(FLEUR_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../")

add_executable(inpgen2
dropDefaults.c
inpgen_help.f90
inpgen.f90
bravais_symm.f90
make_crystal.f90
make_spacegroup.f90
check_mt_radii.f90
make_atom_groups.f90
make_defaults.f90
process_lattice_namelist.f90
atompar.F90
film_sym.f90
make_atomic_defaults.f90
read_inpgen_input.f90
super_check.f90
make_kpoints.f90
read_old_inp.f90
ss_sym.f
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

${FLEUR_SRC}/init/compile_descr.F90
${FLEUR_SRC}/global/sort.f90
${FLEUR_SRC}/math/inv3.f90
${FLEUR_SRC}/io/w_inpXML.f90
${FLEUR_SRC}/kpoints/divi.f
${FLEUR_SRC}/kpoints/bravais.f
${FLEUR_SRC}/kpoints/brzone2.f90
${FLEUR_SRC}/kpoints/kpttet.f
${FLEUR_SRC}/kpoints/kptmop.f
${FLEUR_SRC}/kpoints/ordstar.f
${FLEUR_SRC}/kpoints/fulstar.f
${FLEUR_SRC}/kpoints/tetcon.f
${FLEUR_SRC}/kpoints/kvecon.f
${FLEUR_SRC}/kpoints/kprep.f
${FLEUR_SRC}/kpoints/kptgen_hybrid.f

${FLEUR_SRC}/math/util.F
${FLEUR_SRC}/math/matmul.f
${FLEUR_SRC}/init/strgn_dim.F
${FLEUR_SRC}/init/spgrot.f
#${FLEUR_SRC}/init/convn_dim.f

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
${FLEUR_SRC}/types/types_potden.F90
${FLEUR_SRC}/types/types_enpara.F90
${FLEUR_SRC}/types/types_setup.F90
${FLEUR_SRC}/types/types_usdus.F90
${FLEUR_SRC}/types/types_cdnval.f90
${FLEUR_SRC}/types/types_regionCharges.f90
${FLEUR_SRC}/types/types_dos.f90
${FLEUR_SRC}/types/types_denCoeffsOffdiag.f90
${FLEUR_SRC}/types/types_gpumat.F90

${FLEUR_SRC}/io/nocoInputCheck.F90
${FLEUR_SRC}/eigen/orthoglo.F90
${FLEUR_SRC}/math/ylm4.F90
${FLEUR_SRC}/cdn/slab_dim.f90
${FLEUR_SRC}/cdn/slabgeom.f90
${FLEUR_SRC}/dos/nstm3.f90
${FLEUR_SRC}/cdn/int_21.f90
${FLEUR_SRC}/cdn/int_21lo.f90
${FLEUR_SRC}/cdn_mt/rhomt21.f90
${FLEUR_SRC}/cdn_mt/rhonmt21.f90
${FLEUR_SRC}/xc-pot/gaunt.f90
${FLEUR_SRC}/math/grule.f90
${FLEUR_SRC}/init/ifft235.f
${FLEUR_SRC}/kpoints/gkptwgt.f90
${FLEUR_SRC}/init/local_sym.f
${FLEUR_SRC}/init/lhcal.f
${FLEUR_SRC}/init/ptsym.f
${FLEUR_SRC}/math/gaussp.f
${FLEUR_SRC}/init/gtest.f
${FLEUR_SRC}/init/prp_xcfft_box.f
${FLEUR_SRC}/io/xmlOutput.F90
${FLEUR_SRC}/global/utility.F90
${FLEUR_SRC}/global/radsra.f
${FLEUR_SRC}/global/differ.f
${FLEUR_SRC}/global/soc_sym.f
${FLEUR_SRC}/math/inwint.f
${FLEUR_SRC}/math/outint.f
${FLEUR_SRC}/math/intgr.F90
${FLEUR_SRC}/math/d_wigner.F90
${FLEUR_SRC}/global/find_enpara.f90
${FLEUR_SRC}/io/xsf_io.f90
${FLEUR_SRC}/init/boxdim.f
)

#Set module directories
include_directories("${CMAKE_CURRENT_BINARY_DIR}/modules/inpgen")
set_target_properties(inpgen2 PROPERTIES Fortran_MODULE_DIRECTORY modules/inpgen)
target_compile_definitions(inpgen2 PUBLIC ${FLEUR_DEFINITIONS})
target_link_libraries(inpgen2 ${FLEUR_LIBRARIES})
target_link_libraries(inpgen2 ${FLEUR_LINK_LIBRARIES})
target_link_libraries(inpgen2 juDFT)
target_link_libraries(inpgen2 fleurinput)
