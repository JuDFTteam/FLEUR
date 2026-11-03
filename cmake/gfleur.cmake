#Targets for the Green's-function embedding code (gfleur).
#
#The gfleur executables are built from the full FLEUR source list plus the
#sources in src/gfleur-embedding, replacing only the main program. All
#persistent I/O of gfleur is HDF5-based, so the targets are skipped entirely
#if HDF5 is not available.
#
#As for fleur, there is one target per executable; the definitions and linked
#libraries follow the fleur target in Files_and_Targets.txt.

function(gfleur_configure_target target)
   target_compile_definitions(${target} PUBLIC ${FLEUR_DEFINITIONS} ${FLEUR_PARALLEL_DEFINITIONS})
   target_compile_options(${target} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_PRECISION_OPTION}>)
   target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${FLEUR_COMPILE_OPTIONS}>)
   target_compile_options(${target} PRIVATE -I${CMAKE_BINARY_DIR}/modules/${target})
   set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/${target})
   target_link_libraries(${target} fleurinput juDFT)
   if (FLEUR_COMPILE_EDSOLVER)
      target_link_libraries(${target} EDsolver)
   endif()
   if (FLEUR_COMPILE_MAGMA)
      target_link_libraries(${target} magma)
   endif()
   if (FLEUR_COMPILE_WANNIER90)
      target_link_libraries(${target} Wannier90_lib)
   endif()
   if (FLEUR_USE_MPI)
      if (FLEUR_COMPILE_CHASE)
         target_link_libraries(${target} chase_f)
         set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
      endif()
      if (FLEUR_COMPILE_ELSI)
         target_link_libraries(${target} elsi)
         set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)
      endif()
      if (FLEUR_USE_INTERNAL_ELPA)
         add_dependencies(${target} elpa)
         target_link_libraries(${target} elpa "-lstdc++")
      endif()
      if (FLEUR_COMPILE_SCALAPACK)
         target_link_libraries(${target} scalapack)
      endif()
   endif()
   if (FLEUR_USE_LIBXC AND FLEUR_LIBXC_TARGETS)
      target_link_libraries(${target} ${FLEUR_LIBXC_TARGETS})
   endif()
   target_link_libraries(${target} ${FLEUR_LIBRARIES})
   target_link_libraries(${target} ${FLEUR_LINK_LIBRARIES})
endfunction()

if (FLEUR_USE_HDF5)
   include(src/gfleur-embedding/CMakeLists.txt)
   set(gfleur_SRC ${fleur_SRC} ${gfleur_F90})
   list(REMOVE_ITEM gfleur_SRC src/fleur/main/fleur_prog.F90)
   set_source_files_properties(${gfleur_F90} ${gfprep_F90} PROPERTIES Fortran_FORMAT FREE)
   set(gfprep_SRC ${gfleur_SRC} ${gfprep_F90})
   list(REMOVE_ITEM gfprep_SRC src/gfleur-embedding/main/gfleur.F90)

   #Named 'gfleur_MPI' for a parallel build, like fleur
   add_executable(gfleur ${gfleur_SRC} ${c_filesFleur})
   gfleur_configure_target(gfleur)
   if (FLEUR_USE_MPI)
      set_target_properties(gfleur PROPERTIES RUNTIME_OUTPUT_NAME gfleur_MPI)
   endif()

   #gfprep runs on a single rank
   add_executable(gfprep ${gfprep_SRC} ${c_filesFleur})
   gfleur_configure_target(gfprep)
else()
   message(STATUS "gfleur targets skipped: HDF5 is required")
endif()
