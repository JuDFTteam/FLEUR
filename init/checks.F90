!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_checks
  USE m_juDFT
  IMPLICIT NONE
  private
  public :: check_command_line,check_input_switches
  CONTAINS
    SUBROUTINE check_command_line()
      !Here we check is command line arguments are OK
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      INTEGER:: isize,ierr,irank
#endif
      IF (TRIM(juDFT_string_for_argument("-eig"))=="hdf") THEN
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
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)
      
      IF (irank.EQ.0) THEN
         IF (isize>1) THEN
            IF (TRIM(juDFT_string_for_argument("-eig"))=="-mem") CALL judft_error(&
                "-eig mem cannot be used in parallel mode for Eigenvector IO",hint="Use -eig mpi or -eig hdf instead")
            IF (TRIM(juDFT_string_for_argument("-eig"))=="-da") CALL judft_error(&
                "-eig da cannot be used in parallel mode for Eigenvector IO",hint="Use -eig mpi or -eig hdf instead")
         END IF
      END IF
#endif
    END SUBROUTINE check_command_line

    SUBROUTINE check_input_switches(banddos,vacuum,noco,atoms,input)
      USE m_nocoInputCheck
      USE m_types_fleurinput
      type(t_banddos),INTENT(IN)::banddos
      type(t_vacuum),INTENT(IN) ::vacuum
      type(t_noco),INTENT(IN)   ::noco
      type(t_atoms),INTENT(IN)  ::atoms
      type(t_input),INTENT(IN)  ::input

      integer :: i

           ! Check DOS related stuff (from inped)

     IF ((banddos%ndir.LT.0).AND..NOT.banddos%dos) THEN
        CALL juDFT_error('STOP banddos: the inbuild dos-program  <0'//&
             ' can only be used if dos = true',calledby ="postprocessInput")
     END IF

     IF ((banddos%ndir.LT.0).AND.banddos%dos) THEN
        IF (banddos%e1_dos-banddos%e2_dos.LT.1e-3) THEN
           CALL juDFT_error("STOP banddos: no valid energy window for "//&
                "internal dos-program",calledby ="postprocessInput")
        END IF
        IF (banddos%sig_dos.LT.0) THEN
           CALL juDFT_error("STOP DOS: no valid broadening (sig_dos) for "//&
                "internal dos-PROGRAM",calledby ="postprocessInput")
        END IF
     END IF

     IF (banddos%vacdos) THEN
        IF (.NOT.banddos%dos) THEN
           CALL juDFT_error("STOP DOS: only set vacdos = .true. if dos = .true.",calledby ="postprocessInput")
        END IF
        IF (.NOT.vacuum%starcoeff.AND.(vacuum%nstars.NE.1))THEN
           CALL juDFT_error("STOP banddos: if stars = f set vacuum=1",calledby ="postprocessInput")
        END IF
        IF (vacuum%layers.LT.1) THEN
           CALL juDFT_error("STOP DOS: specify layers if vacdos = true",calledby ="postprocessInput")
        END IF
        DO i=1,vacuum%layers
           IF (vacuum%izlay(i,1).LT.1) THEN
              CALL juDFT_error("STOP DOS: all layers must be at z>0",calledby ="postprocessInput")
           END IF
        END DO
     END IF
     IF (noco%l_noco) CALL nocoInputCheck(atoms,input,vacuum,noco)
   END SUBROUTINE check_input_switches

  END MODULE m_checks
