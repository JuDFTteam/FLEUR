!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_checks
  USE m_juDFT
  CONTAINS
    SUBROUTINE check_command_line()
      !Here we check is command line arguments are OK
      IMPLICIT NONE
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER:: isize,ierr,irank
#endif
      IF (juDFT_was_argument("-hdf")) THEN
#ifndef CPP_HDF
         CALL judft_error("HDF5 cannot be used for Eigenvector IO",&
              hint="You compiled without support for HDF5. Please use another mode")
#endif
#ifdef CPP_MPI
#ifndef CPP_HDFMPI
         CALL judft_error("HDF5 cannot be used in parallel mode for Eigenvector IO",&
              hint="Your HDF5 library does not support parallel IO" )
#endif
#endif
      ENDIF
      !Check for IO options not available in parallel
#ifdef CPP_MPI
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,irank,ierr)
      IF (isize>1) THEN
         IF (juDFT_was_argument("-mem")) CALL judft_error(&
              "-mem cannot be used in parallel mode for Eigenvector IO",hint="Use -mpi or -hdf instead")
         IF (juDFT_was_argument("-da")) CALL judft_error(&
              "-da cannot be used in parallel mode for Eigenvector IO",hint="Use -mpi or -hdf instead")
      END IF
#endif
    END SUBROUTINE check_command_line
  END MODULE m_checks
