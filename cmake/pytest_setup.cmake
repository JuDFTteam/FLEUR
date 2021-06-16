#this file stores settings to be used in the testing-system

#some test need specific FLEUR features
if (NOT FLEUR_USE_HDF5)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} hdf")
endif()
if (NOT FLEUR_USE_LIBXC)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} libxc")
endif()
if (NOT FLEUR_USE_FFTMKL)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} fftmkl")
endif()
if (NOT FLEUR_USE_FFTW)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} fftw")
endif()
if (NOT FLEUR_USE_SPFFT)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} spfft")
endif()
if (NOT FLEUR_USE_WANN)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} wannier")
endif()
if (NOT FLEUR_USE_WANN4)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} wannier")
endif()
if (NOT FLEUR_USE_WANN5)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} wannier")
endif()
if (NOT FLEUR_USE_MAGMA)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} magma")
endif()
if (NOT FLEUR_USE_EDSOLVER)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} edsolver")
endif()
if (NOT FLEUR_USE_CUSOLVER)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} cusolver")
endif()
if (NOT FLEUR_USE_PROG_THREAD)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} progthread")
endif()
if (NOT FLEUR_USE_ELPA)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} elpa")
endif()
if (NOT FLEUR_USE_ELPA_ONENODE)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} elpaonenode")
endif()
if (NOT FLEUR_USE_CHASE)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} chase")
endif()
if (NOT FLEUR_USE_SCALAPACK)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} scalapack")
endif()
if (NOT FLEUR_USE_GPU)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} gpu")
endif()
if (NOT FLEUR_USE_MPI)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} mpionly")
endif()
if (NOT FLEUR_USE_OPENMP)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} openmponly")
endif()
#write file
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/pytest_incl.py CONTENT "sourcedir=${CMAKE_SOURCE_DIR}\nbuilddir=${CMAKE_BINARY_DIR}\nexcl_flags=\"${PYTEST_TEST_EXCL_FLAGS}\"\n")

get_filename_component(BUILD_DIR ${CMAKE_BINARY_DIR} NAME)
#write build script
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/run_tests.sh CONTENT
"#!/usr/bin/env bash
PYTEST_ADDOPTS=\"../tests/new_pytest_system --build_dir=../../${BUILD_DIR} $PYTEST_ADD_OPTS\"
if [[ -z \"\${juDFT_PYTHON}\" ]]; then
  PYTEST_ADDOPTS=$PYTEST_ADDOPTS pytest \"$@\"
else
  PYTEST_ADDOPTS=$PYTEST_ADDOPTS $juDFT_PYTHON -m pytest \"$@\"
fi")
add_custom_target(pytest ALL
                  COMMAND chmod +x run_tests.sh
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  COMMENT "Making test script executable")


