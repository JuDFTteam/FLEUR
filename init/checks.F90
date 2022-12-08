!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_checks
  USE m_juDFT
  USE m_types
  IMPLICIT NONE
  private
  public :: check_command_line,check_input_switches
  CONTAINS
    SUBROUTINE check_command_line(fmpi)
      !Here we check is command line arguments are OK
#ifdef CPP_MPI
      USE mpi
      INTEGER:: isize,ierr,irank
#endif
      TYPE(t_mpi), INTENT(INOUT):: fmpi
      IF (TRIM(juDFT_string_for_argument("-eig"))=="hdf") THEN
#ifndef CPP_HDF
         CALL judft_error("HDF5 cannot be used for Eigenvector IO",&
              hint="You compiled without support for HDF5. Please use another mode")
#endif
#ifdef CPP_MPI
#ifndef CPP_HDFMPI
         IF (fmpi%irank.GT.1) THEN
            CALL judft_error("HDF5 cannot be used in parallel mode for Eigenvector IO",&
                 hint="Your HDF5 library does not support parallel IO" )
         END IF
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

    SUBROUTINE check_input_switches(banddos,vacuum,noco,atoms,input,sym,kpts)
      USE m_nocoInputCheck
      USE m_types_fleurinput
      USE m_constants
      type(t_banddos),INTENT(IN)::banddos
      type(t_vacuum),INTENT(IN) ::vacuum
      type(t_noco),INTENT(IN)   ::noco
      type(t_atoms),INTENT(IN)  ::atoms
      type(t_input),INTENT(IN)  ::input
      type(t_sym),INTENT(IN)    :: sym
      type(t_kpts),INTENT(IN)   :: kpts

      integer :: i,n,na
      real :: maxpos,minpos

     ! Check DOS related stuff (from inped)
     IF(banddos%l_jDOS.AND..NOT.noco%l_noco) THEN
        CALL juDFT_error("jDOS+collinear is not implemented at the moment.",&
                         hint="If you need this feature please create an issue on the fleur git")
     ENDIF

      IF(atoms%n_opc>0.AND..NOT.noco%l_soc) THEN
         CALL juDFT_error("Orbital polarization corrections without spin-orbit coupling have no effect",&
                          calledby="check_input_switches")
      ENDIF

     IF (banddos%vacdos) THEN
        IF (.NOT.banddos%dos) THEN
           CALL juDFT_error("STOP DOS: only set vacdos = .true. if dos = .true.",calledby ="check_input_switches")
        END IF
        IF (.NOT.banddos%starcoeff.AND.(banddos%nstars.NE.1))THEN
           CALL juDFT_error("STOP banddos: if stars = f set vacuum=1",calledby ="check_input_switches")
        END IF
        IF (banddos%layers.LT.1) THEN
           CALL juDFT_error("STOP DOS: specify layers if vacdos = true",calledby ="check_input_switches")
        END IF
        DO i=1,banddos%layers
           IF (banddos%izlay(i,1).LT.1) THEN
              CALL juDFT_error("STOP DOS: all layers must be at z>0",calledby ="check_input_switches")
           END IF
        END DO
     END IF

     IF((input%gw.EQ.2).AND.(kpts%kptsKind.NE.KPTS_KIND_SPEX_MESH)) THEN
        CALL juDFT_warn('Chosen k-point set is not compatible to this GW step.', calledby='check_input_switches')
     END IF

     IF (noco%l_noco) CALL nocoInputCheck(atoms,input,sym,vacuum,noco)

     !In film case check centering of film
     if ( input%film ) then
        IF ((input%f_level.GT.0.).AND.input%l_f) THEN
           call judft_warn("Enhanced forces are not implemented for film calculations.",hint="Set the f_level tag to 0.")
        END IF
       maxpos=0.0;minpos=0.0
       DO n=1,atoms%ntype
         na=sum(atoms%neq(:n-1))
         maxpos=max(maxpos,maxval(atoms%pos(3,na+1:na+atoms%neq(n)))+atoms%rmt(n))
         minpos=max(minpos,maxval(-1.*atoms%pos(3,na+1:na+atoms%neq(n)))+atoms%rmt(n))
       ENDDO
       if (abs(maxpos-minpos)>2.0) call judft_warn("Your film setup is not centered around zero",hint="Using a non-centred setup can lead to numerical problems. Please check your setup an try to ensure that the center of your film is at z=0")
     endif

     IF(banddos%unfoldBand.AND.noco%l_soc.AND..NOT.noco%l_noco) THEN
        IF(banddos%unfoldUseOlap) THEN
           CALL juDFT_error("Band unfolding for 2nd variation SOC is only implemented without incorporating the overlap matrix.", hint="You have to set the optional /output/unfoldingBand/@useOlap switch to F.", calledby ="check_input_switches")
        END IF
     END IF

   END SUBROUTINE check_input_switches

  END MODULE m_checks
