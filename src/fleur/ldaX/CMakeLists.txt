set(fleur_F77 ${fleur_F77}
)
set(fleur_F90 ${fleur_F90}
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/sgaunt.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/u_setup.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/uj2f.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/umtx.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/doubleCounting.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/dftUPotential.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/coulombPotential.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/v_ham.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/add_selfen.F90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/crystalfield.F90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/hubbard1_setup.F90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/hubbard1Distance.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/opc_setup.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/slater.f90
${CMAKE_SOURCE_DIR}/src/fleur/ldaX/types_selfen.F90
)
