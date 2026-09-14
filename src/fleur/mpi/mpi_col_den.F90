!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_col_den
  !
  ! collect all data calculated in cdnval on different pe's on pe 0
  !
  ! for some data also spread them back onto all pe's (Jan. 2019  U.Alekseeva)
  !
#ifdef CPP_MPI
   use mpi
#endif
   implicit none
CONTAINS
  SUBROUTINE mpi_col_den(fmpi,sphhar,atoms ,stars,vacuum,input,noco,jspin,dos,vacdos,&
                         results,den,mcd,slab,orbcomp,jDOS)

    USE m_types
    USE m_constants
    USE m_juDFT
    use m_types_mcd
    use m_types_slab
    use m_types_orbcomp
    use m_types_jDOS
    use m_types_vacdos
    IMPLICIT NONE

    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_mpi),INTENT(IN)       :: fmpi

    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_potden),INTENT(INOUT) :: den
    ! ..
    ! ..  Scalar Arguments ..
    INTEGER, INTENT (IN) :: jspin
    ! ..
    ! ..  Array Arguments ..

    TYPE (t_dos),               INTENT(INOUT) :: dos
    TYPE (t_vacdos),            INTENT(INOUT) :: vacdos
    !TYPE (t_regionCharges), OPTIONAL, INTENT(INOUT) :: regCharges
    TYPE (t_mcd),       OPTIONAL, INTENT(INOUT) :: mcd
    TYPE (t_slab),      OPTIONAL, INTENT(INOUT) :: slab
    TYPE (t_orbcomp),   OPTIONAL, INTENT(INOUT) :: orbcomp
    TYPE (t_jDOS),      OPTIONAL, INTENT(INOUT) :: jDOS
    ! ..
    ! ..  Local Scalars ..
    INTEGER :: n, i
    ! ..
    ! ..  Local Arrays ..
    INTEGER :: ierr
    COMPLEX, ALLOCATABLE :: c_b(:)
    REAL,    ALLOCATABLE :: r_b(:)
    INTEGER, ALLOCATABLE :: i_b(:)
    ! ..
    ! ..  External Subroutines
#ifdef CPP_MPI
    CALL timestart("mpi_col_den")

    ! -> Collect den%pw(:,jspin)
    n = stars%ng3
    ALLOCATE(c_b(n))
    CALL MPI_ALLREDUCE(den%pw(:,jspin),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL zcopy(n, c_b, 1, den%pw(:,jspin), 1)
    DEALLOCATE (c_b)

    IF (input%film) THEN
       ! -> Collect den%vac(:,:,:,jspin)
       n=size(den%vac(:,:,:,jspin))
       ALLOCATE(c_b(n))
       CALL MPI_REDUCE(den%vac(:,:,:,jspin),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) CALL zcopy(n, c_b, 1, den%vac(:,:,:,jspin), 1)
       DEALLOCATE (c_b)
    ENDIF

    
    !collect DOS stuff
    if (dos%l_initialized) then
      n = SIZE(dos%jsym,1)*SIZE(dos%jsym,2)
      ALLOCATE(i_b(n))
      CALL MPI_REDUCE(dos%jsym(:,:,jspin),i_b,n,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) THEN
         DO i = 1, SIZE(dos%jsym,2)
            dos%jsym(:,i,jspin) = i_b((i-1)*SIZE(dos%jsym,1)+1:i*SIZE(dos%jsym,1))
         END DO
      END IF
      DEALLOCATE (i_b)

      n = SIZE(dos%qis,1)*SIZE(dos%qis,2)
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(dos%qis(:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, dos%qis(:,:,jspin), 1)
      DEALLOCATE (r_b)

      n = SIZE(dos%qTot,1)*SIZE(dos%qTot,2)
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(dos%qTot(:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, dos%qTot(:,:,jspin), 1)
      DEALLOCATE (r_b)

      n = SIZE(dos%qal,1)*SIZE(dos%qal,2)*SIZE(dos%qal,3)*SIZE(dos%qal,4)
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(dos%qal(0:,:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, dos%qal(0:,:,:,:,jspin), 1)
      DEALLOCATE (r_b)
    END IF
    if (vacdos%l_initialized) then
      
      n = SIZE(vacdos%qvac,1)*SIZE(vacdos%qvac,2)*SIZE(vacdos%qvac,3)
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(vacdos%qvac(:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, vacdos%qvac(:,:,:,jspin), 1)
      DEALLOCATE (r_b)

      n = SIZE(vacdos%qvlay,1)*SIZE(vacdos%qvlay,2)*SIZE(vacdos%qvlay,3)*SIZE(vacdos%qvlay,4)
      ALLOCATE(r_b(n))
      CALL MPI_REDUCE(vacdos%qvlay(:,:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, vacdos%qvlay(:,:,:,:,jspin), 1)
      DEALLOCATE (r_b)

      n = SIZE(vacdos%qstars,1)*SIZE(vacdos%qstars,2)*SIZE(vacdos%qstars,3)*SIZE(vacdos%qstars,4)*SIZE(vacdos%qstars,5)
      ALLOCATE(c_b(n))
      CALL MPI_REDUCE(vacdos%qstars(:,:,:,:,:,jspin),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) CALL zcopy(n, c_b, 1, vacdos%qstars(:,:,:,:,:,jspin), 1)
      DEALLOCATE (c_b)
    endif
    ! Collect mcd%mcd
    IF (PRESENT(mcd)) THEN
      if (mcd%l_initialized) then
       n = SIZE(mcd%mcd,1)*SIZE(mcd%mcd,2)*SIZE(mcd%mcd,3)*SIZE(mcd%mcd,4)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(mcd%mcd(:,:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, mcd%mcd(:,:,:,:,jspin), 1)
       DEALLOCATE (r_b)
      endif
    END IF

    ! Collect slab - qintsl and qmtsl
    IF (PRESENT(slab)) THEN
      if (slab%l_initialized) then
       n = SIZE(slab%qintsl,1)*SIZE(slab%qintsl,2)*SIZE(slab%qintsl,3)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(slab%qintsl(:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, slab%qintsl(:,:,:,jspin), 1)
       DEALLOCATE (r_b)

       n = SIZE(slab%qmtsl,1)*SIZE(slab%qmtsl,2)*SIZE(slab%qmtsl,3)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(slab%qmtsl(:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, slab%qmtsl(:,:,:,jspin), 1)
       DEALLOCATE (r_b)
      endif
    END IF

    ! Collect orbcomp - comp and qmtp
    IF (PRESENT(orbcomp)) THEN
      if (orbcomp%l_initialized) then
       n = SIZE(orbcomp%comp,1)*SIZE(orbcomp%comp,2)*SIZE(orbcomp%comp,3)*SIZE(orbcomp%comp,4)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(orbcomp%comp(:,:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, orbcomp%comp(:,:,:,:,jspin), 1)
       DEALLOCATE (r_b)

       n = SIZE(orbcomp%qmtp,1)*SIZE(orbcomp%qmtp,2)*SIZE(orbcomp%qmtp,3)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(orbcomp%qmtp(:,:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, orbcomp%qmtp(:,:,:,jspin), 1)
       DEALLOCATE (r_b)
      endif 
    END IF

    !+jDOS
    IF(PRESENT(jDOS)) THEN
      IF(jdos%l_initialized.and.jspin.EQ.1) THEN

        n = SIZE(jDOS%comp)
        ALLOCATE(r_b(n))
        CALL MPI_REDUCE(jDOS%comp,r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF(fmpi%irank.EQ.0) CALL dcopy(n,r_b,1,jDOS%comp,1)
        DEALLOCATE(r_b)

        n = SIZE(jDOS%qmtp)
        ALLOCATE(r_b(n))
        CALL MPI_REDUCE(jDOS%qmtp,r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF(fmpi%irank.EQ.0) CALL dcopy(n,r_b,1,jDOS%qmtp,1)
        DEALLOCATE(r_b)

        n = SIZE(jDOS%occ)
        ALLOCATE(r_b(n))
        CALL MPI_REDUCE(jDOS%occ,r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        IF(fmpi%irank.EQ.0) CALL dcopy(n,r_b,1,jDOS%occ,1)
        DEALLOCATE(r_b)

      ENDIF
    ENDIF
    !-jDOS

    ! -> Collect force
    IF (input%l_f) THEN
       n=3*atoms%ntype
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(results%force(1,1,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) CALL dcopy(n, r_b, 1, results%force(1,1,jspin), 1)
       DEALLOCATE (r_b)
    ENDIF

    ! -> Optional the LO-coefficients: aclo,bclo,enerlo,cclo,acnmt,bcnmt,ccnmt
    IF (atoms%nlod.GE.1) THEN

       n=atoms%nlod*atoms%ntype
       ALLOCATE (r_b(n))
#if false       
       IF (PRESENT(regCharges)) THEN
         CALL MPI_ALLREDUCE(regCharges%enerlo(:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
         CALL dcopy(n, r_b, 1, regCharges%enerlo(:,:,jspin), 1)
         CALL MPI_ALLREDUCE(regCharges%sqlo(:,:,jspin),r_b,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
         CALL dcopy(n, r_b, 1, regCharges%sqlo(:,:,jspin), 1)
       END IF
#endif       
       DEALLOCATE (r_b)

       
    ENDIF

    
    ! -> Collect the noco stuff:
    IF ( noco%l_noco .AND. jspin.EQ.1 ) THEN

       n = stars%ng3
       ALLOCATE(c_b(n))
       CALL MPI_REDUCE(den%pw(:,3),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) THEN
          den%pw(:,3)=RESHAPE(c_b,(/n/))
       ENDIF
       DEALLOCATE (c_b)
       !
       IF (input%film) THEN
          n=size(den%vac(:,:,:,3))
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(den%vac(:,:,:,3),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
          IF (fmpi%irank.EQ.0) THEN
             CALL zcopy(n, c_b, 1, den%vac(:,:,:,3), 1)
          ENDIF
          DEALLOCATE (c_b)
       ENDIF ! input%film


     
    ENDIF   ! noco

    !+lda+U
    IF ( atoms%n_u.GT.0 ) THEN
       n = 49*atoms%n_u
       ALLOCATE(c_b(n))
       CALL MPI_REDUCE(den%mmpMat(:,:,1:atoms%n_u,jspin),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (fmpi%irank.EQ.0) THEN
          CALL zcopy(n, c_b, 1, den%mmpMat(:,:,1:atoms%n_u,jspin), 1)
       ENDIF
       DEALLOCATE (c_b)
       IF(noco%l_mperp.AND.jspin.EQ.1) THEN
         n = 49*atoms%n_u
         ALLOCATE(c_b(n))
         CALL MPI_REDUCE(den%mmpMat(:,:,1:atoms%n_u,3),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
         IF (fmpi%irank.EQ.0) THEN
            CALL zcopy(n, c_b, 1, den%mmpMat(:,:,1:atoms%n_u,3), 1)
         ENDIF
         DEALLOCATE (c_b)
       ENDIF
    ENDIF
    !-lda+U
       
    !+LDA+V
    IF ( atoms%n_v.GT.0 ) THEN
       ALLOCATE(c_b(size(den%nIJ_llp_mmp)))
       CALL MPI_REDUCE(den%nIJ_llp_mmp,c_b,size(den%nIJ_llp_mmp),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       if (fmpi%irank.EQ.0) den%nIJ_llp_mmp=reshape(c_b,shape(den%nIJ_llp_mmp))
       DEALLOCATE(c_b)
    END IF
    !-LDA+V

    !+lda+OP
    IF ( atoms%n_opc.GT.0 ) THEN
      n = 49*atoms%n_opc
      ALLOCATE(c_b(n))
      CALL MPI_REDUCE(den%mmpMat(:,:,atoms%n_u+atoms%n_hia+1:,jspin),c_b,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
      IF (fmpi%irank.EQ.0) THEN
         CALL zcopy(n, c_b, 1, den%mmpMat(:,:,atoms%n_u+atoms%n_hia+1:,jspin), 1)
      ENDIF
      DEALLOCATE (c_b)
   ENDIF
   !-lda+U

    CALL timestop("mpi_col_den")

#endif

  END SUBROUTINE mpi_col_den
END MODULE m_mpi_col_den
