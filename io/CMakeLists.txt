set(fleur_F77 ${fleur_F77}
io/cdn_read.F
)
set(fleur_F90 ${fleur_F90}
io/io_matrix.F90
io/io_hybinp.F90
io/relax_io.F90
io/eig66_da.F90
io/eig66_data.F90
io/eig66_hdf.F90
io/eig66_io.F90
io/eig66_mem.F90
io/eig66_mpi.F90
io/nocoInputCheck.F90
io/loddop.f90
io/cdnpot_io_hdf.F90
io/cdnpot_io_common.F90
io/cdn_io.F90
io/pot_io.F90
io/banddos_io.F90
#io/r_inpXML.F90
io/writeBasis.F90
io/writeOutHeader.F90
io/wrtdop.f90
io/w_inpXML.f90
io/writeOutParameters.f90
io/xsf_io.f90
io/xmlOutput.F90
io/hubbard1_io.f90
io/generic_txtio.f90
)
if (FLEUR_USE_HDF5)
set(fleur_F90 ${fleur_F90}
io/iomatrix_hdf.F90)
endif()

