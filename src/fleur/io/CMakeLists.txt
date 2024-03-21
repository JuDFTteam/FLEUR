set(fleur_F77 ${fleur_F77}
src/fleur/io/cdn_read.F
)
set(fleur_F90 ${fleur_F90}
src/fleur/io/io_matrix.F90
src/fleur/io/io_hybrid.F90
src/fleur/io/relax_io.F90
src/fleur/io/eig66_da.F90
src/fleur/io/eig66_data.F90
src/fleur/io/eig66_hdf.F90
src/fleur/io/eig66_io.F90
src/fleur/io/eig66_mem.F90
src/fleur/io/eig66_mpi.F90
src/fleur/io/nocoInputCheck.F90
src/fleur/io/loddop.f90
src/fleur/io/cdnpot_io_hdf.F90
src/fleur/io/cdnpot_io_common.F90
src/fleur/io/cdn_io.F90
src/fleur/io/pot_io.F90
src/fleur/io/banddos_io.F90
src/fleur/io/greensf_io.F90
#src/fleur/io/r_inpXML.F90
src/fleur/io/writeBasis.F90
src/fleur/io/writeOutHeader.F90
src/fleur/io/wrtdop.f90
src/fleur/io/w_inpXML.f90
src/fleur/io/writeOutParameters.f90
src/fleur/io/xsf_io.f90
src/fleur/io/xmlOutput.F90
src/fleur/io/hubbard1_io.f90
src/fleur/io/generic_txtio.f90
src/fleur/io/cfOutput_hdf.F90
src/fleur/io/writeCFOutput.F90
)
if (FLEUR_USE_HDF5)
set(fleur_F90 ${fleur_F90}
src/fleur/io/iomatrix_hdf.F90)
endif()
