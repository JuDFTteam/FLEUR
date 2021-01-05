#this file stores settings to be used in the testing-system

#some test need specific FLEUR features
if (NOT FLEUR_USE_HDF5)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} hdf5")
endif()
if (NOT FLEUR_USE_LIBXC)
    set(PYTEST_TEST_EXCL_FLAGS "${PYTEST_TEST_EXCL_FLAGS} libxc")
endif()

#write file
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/pytest_incl.py CONTENT "sourcedir=${CMAKE_SOURCE_DIR}\nbuilddir=${CMAKE_BINARY_DIR}\nexcl_flags=\"${PYTEST_TEST_EXCL_FLAGS}\"\n")
