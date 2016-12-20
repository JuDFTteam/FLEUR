set(fleur_F77 ${fleur_F77}
)
if (${FLEUR_USE_HDF5})
set(fleur_F90 ${fleur_F90}
juDFT/hdf_accessprp.F90
juDFT/hdf_tools.F90
juDFT/hdf_tools_attr.F90
juDFT/hdf_tools_exists.F90
juDFT/hdf_tools_misc.F90
juDFT/hdf_tools_rw.F90
juDFT/hdf_tools_rw_var.F90
juDFT/hdf_tools_stride.F90)
endif()
set(fleur_F90 ${fleur_F90}
juDFT/init.F90
juDFT/juDFT.F90
juDFT/stop.F90
juDFT/time.F90
juDFT/args.F90
juDFT/sysinfo.F90
)
