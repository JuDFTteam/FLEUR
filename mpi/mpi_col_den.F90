!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_col_den
  !
  ! collect all data calculated in cdnval on different pe's on pe 0
  !
CONTAINS
  SUBROUTINE mpi_col_den(mpi,sphhar,atoms,oneD,stars,vacuum,input,noco,jspin,regCharges,dos,&
                         results,denCoeffs,orb,denCoeffsOffdiag,den,n_mmp,mcd,slab,orbcomp)

#include"cpp_double.h"
    USE m_types
    USE m_constants
    USE m_juDFT
    IMPLICIT NONE

    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_potden),INTENT(INOUT) :: den
    INCLUDE 'mpif.h'
    ! ..
    ! ..  Scalar Arguments ..
    INTEGER, INTENT (IN) :: jspin
    ! ..
    ! ..  Array Arguments ..
    COMPLEX,INTENT(INOUT) :: n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u)
    TYPE (t_orb),               INTENT(INOUT) :: orb
    TYPE (t_denCoeffs),         INTENT(INOUT) :: denCoeffs
    TYPE (t_denCoeffsOffdiag),  INTENT(INOUT) :: denCoeffsOffdiag
    TYPE (t_regionCharges),     INTENT(INOUT) :: regCharges
    TYPE (t_dos),               INTENT(INOUT) :: dos
    TYPE (t_mcd),     OPTIONAL, INTENT(INOUT) :: mcd
    TYPE (t_slab),    OPTIONAL, INTENT(INOUT) :: slab
    TYPE (t_orbcomp), OPTIONAL, INTENT(INOUT) :: orbcomp
    ! ..
    ! ..  Local Scalars ..
    INTEGER :: n, i
    ! ..
    ! ..  Local Arrays ..
    INTEGER :: ierr(3)
    COMPLEX, ALLOCATABLE :: c_b(:)
    REAL,    ALLOCATABLE :: r_b(:)
    INTEGER, ALLOCATABLE :: i_b(:)
    ! ..
    ! ..  External Subroutines
    EXTERNAL CPP_BLAS_scopy,CPP_BLAS_ccopy,MPI_REDUCE

    CALL timestart("mpi_col_den")

    ! -> Collect den%pw(:,jspin)
    n = stars%ng3
    ALLOCATE(c_b(n))
    CALL MPI_REDUCE(den%pw(:,jspin),c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_ccopy(n, c_b, 1, den%pw(:,jspin), 1)
    DEALLOCATE (c_b)

    ! -> Collect den%vacxy(:,:,:,jspin)
    IF (input%film) THEN
       n = vacuum%nmzxy*(oneD%odi%n2d-1)*2
       ALLOCATE(c_b(n))
       CALL MPI_REDUCE(den%vacxy(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_ccopy(n, c_b, 1, den%vacxy(:,:,:,jspin), 1)
       DEALLOCATE (c_b)

       ! -> Collect den%vacz(:,:,jspin)
       n = vacuum%nmz*2
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(den%vacz(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, den%vacz(:,:,jspin), 1)
       DEALLOCATE (r_b)
    ENDIF

    ! -> Collect uu(),ud() and dd()
    n = (atoms%lmaxd+1)*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(denCoeffs%uu(0:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%uu(0:,:,jspin), 1)
    CALL MPI_REDUCE(denCoeffs%du(0:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%du(0:,:,jspin), 1)
    CALL MPI_REDUCE(denCoeffs%dd(0:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%dd(0:,:,jspin), 1)
    DEALLOCATE (r_b)

    !--> Collect uunmt,udnmt,dunmt,ddnmt
    n = (((atoms%lmaxd*(atoms%lmaxd+3))/2)+1)*sphhar%nlhd*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(denCoeffs%uunmt(0:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%uunmt(0:,:,:,jspin), 1)
    CALL MPI_REDUCE(denCoeffs%udnmt(0:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%udnmt(0:,:,:,jspin), 1)
    CALL MPI_REDUCE(denCoeffs%dunmt(0:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%dunmt(0:,:,:,jspin), 1)
    CALL MPI_REDUCE(denCoeffs%ddnmt(0:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%ddnmt(0:,:,:,jspin), 1)
    DEALLOCATE (r_b)

    !--> ener & sqal
    n=4*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(regCharges%ener(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, regCharges%ener(:,:,jspin), 1)
    CALL MPI_REDUCE(regCharges%sqal(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, regCharges%sqal(:,:,jspin), 1)
    DEALLOCATE (r_b)

    !--> svac & pvac
    IF ( input%film ) THEN
       n=2
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(regCharges%svac(:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, regCharges%svac(:,jspin), 1)
       CALL MPI_REDUCE(regCharges%pvac(:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, regCharges%pvac(:,jspin), 1)
       DEALLOCATE (r_b)
    ENDIF

    !collect DOS stuff
    n = SIZE(dos%jsym,1)*SIZE(dos%jsym,2)
    ALLOCATE(i_b(n))
    CALL MPI_REDUCE(dos%jsym(:,:,jspin),i_b,n,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       DO i = 1, SIZE(dos%jsym,2)
          dos%jsym(:,i,jspin) = i_b((i-1)*SIZE(dos%jsym,1)+1:i*SIZE(dos%jsym,1))
       END DO
    END IF
    DEALLOCATE (i_b)

    n = SIZE(dos%ksym,1)*SIZE(dos%ksym,2)
    ALLOCATE(i_b(n))
    CALL MPI_REDUCE(dos%ksym(:,:,jspin),i_b,n,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       DO i = 1, SIZE(dos%ksym,2)
          dos%ksym(:,i,jspin) = i_b((i-1)*SIZE(dos%ksym,1)+1:i*SIZE(dos%ksym,1))
       END DO
    END IF
    DEALLOCATE (i_b)

    n = SIZE(dos%qis,1)*SIZE(dos%qis,2)
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(dos%qis(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, dos%qis(:,:,jspin), 1)
    DEALLOCATE (r_b)

    n = SIZE(dos%qal,1)*SIZE(dos%qal,2)*SIZE(dos%qal,3)*SIZE(dos%qal,4)
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(dos%qal(0:,:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, dos%qal(0:,:,:,:,jspin), 1)
    DEALLOCATE (r_b)

    n = SIZE(dos%qvac,1)*SIZE(dos%qvac,2)*SIZE(dos%qvac,3)
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(dos%qvac(:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, dos%qvac(:,:,:,jspin), 1)
    DEALLOCATE (r_b)

    n = SIZE(dos%qvlay,1)*SIZE(dos%qvlay,2)*SIZE(dos%qvlay,3)*SIZE(dos%qvlay,4)
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(dos%qvlay(:,:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, dos%qvlay(:,:,:,:,jspin), 1)
    DEALLOCATE (r_b)

    n = SIZE(dos%qstars,1)*SIZE(dos%qstars,2)*SIZE(dos%qstars,3)*SIZE(dos%qstars,4)*SIZE(dos%qstars,5)
    ALLOCATE(c_b(n))
    CALL MPI_REDUCE(dos%qstars(:,:,:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) CALL CPP_BLAS_ccopy(n, c_b, 1, dos%qstars(:,:,:,:,:,jspin), 1)
    DEALLOCATE (c_b)

    ! Collect mcd%mcd
    IF (PRESENT(mcd)) THEN
       n = SIZE(mcd%mcd,1)*SIZE(mcd%mcd,2)*SIZE(mcd%mcd,3)*SIZE(mcd%mcd,4)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(mcd%mcd(:,:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, mcd%mcd(:,:,:,:,jspin), 1)
       DEALLOCATE (r_b)
    END IF

    ! Collect slab - qintsl and qmtsl
    IF (PRESENT(slab)) THEN
       n = SIZE(slab%qintsl,1)*SIZE(slab%qintsl,2)*SIZE(slab%qintsl,3)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(slab%qintsl(:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, slab%qintsl(:,:,:,jspin), 1)
       DEALLOCATE (r_b)

       n = SIZE(slab%qmtsl,1)*SIZE(slab%qmtsl,2)*SIZE(slab%qmtsl,3)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(slab%qmtsl(:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, slab%qmtsl(:,:,:,jspin), 1)
       DEALLOCATE (r_b)
    END IF

    ! Collect orbcomp - comp and qmtp
    IF (PRESENT(orbcomp)) THEN
       n = SIZE(orbcomp%comp,1)*SIZE(orbcomp%comp,2)*SIZE(orbcomp%comp,3)*SIZE(orbcomp%comp,4)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(orbcomp%comp(:,:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, orbcomp%comp(:,:,:,:,jspin), 1)
       DEALLOCATE (r_b)

       n = SIZE(orbcomp%qmtp,1)*SIZE(orbcomp%qmtp,2)*SIZE(orbcomp%qmtp,3)
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(orbcomp%qmtp(:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, orbcomp%qmtp(:,:,:,jspin), 1)
       DEALLOCATE (r_b)
    END IF

    ! -> Collect force
    IF (input%l_f) THEN
       n=3*atoms%ntype
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(results%force(1,1,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) CALL CPP_BLAS_scopy(n, r_b, 1, results%force(1,1,jspin), 1)
       DEALLOCATE (r_b)
    ENDIF

    ! -> Optional the LO-coefficients: aclo,bclo,enerlo,cclo,acnmt,bcnmt,ccnmt
    IF (atoms%nlod.GE.1) THEN

       n=atoms%nlod*atoms%ntype 
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(denCoeffs%aclo(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%aclo(:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(denCoeffs%bclo(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%bclo(:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(regCharges%enerlo(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, regCharges%enerlo(:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(regCharges%sqlo(:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, regCharges%sqlo(:,:,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

       n = atoms%nlod * atoms%nlod * atoms%ntype
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(denCoeffs%cclo(:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%cclo(:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

       n = (atoms%lmaxd+1) * atoms%ntype * atoms%nlod * sphhar%nlhd
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(denCoeffs%acnmt(0:,:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%acnmt(0:,:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(denCoeffs%bcnmt(0:,:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%bcnmt(0:,:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

       n = atoms%ntype * sphhar%nlhd * atoms%nlod**2
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(denCoeffs%ccnmt(:,:,:,:,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, denCoeffs%ccnmt(:,:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

    ENDIF

    ! ->  Now the SOC - stuff: orb, orblo and orblo
    IF (noco%l_soc) THEN
       ! orb
       n=(atoms%lmaxd+1)*(2*atoms%lmaxd+1)*atoms%ntype
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(orb%uu(:,:,:,jspin),r_b,n,CPP_MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orb%uu(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%dd(:,:,:,jspin),r_b,n,CPP_MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orb%dd(:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

       ALLOCATE (c_b(n))
       CALL MPI_REDUCE(orb%uup(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%uup(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%ddp(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%ddp(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%uum(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%uum(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%ddm(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%ddm(:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (c_b)

       n = atoms%nlod * (2*atoms%llod+1) * atoms%ntype
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(orb%uulo(:,:,:,jspin),r_b,n,CPP_MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orb%uulo(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%dulo(:,:,:,jspin),r_b,n,CPP_MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orb%dulo(:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

       ALLOCATE (c_b(n))
       CALL MPI_REDUCE(orb%uulop(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%uulop(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%dulop(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%dulop(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%uulom(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%uulom(:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%dulom(:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%dulom(:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (c_b)

       n = atoms%nlod * atoms%nlod * (2*atoms%llod+1) * atoms%ntype
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(orb%z(:,:,:,:,jspin),r_b,n,CPP_MPI_REAL, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, orb%z(:,:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

       ALLOCATE (c_b(n))
       CALL MPI_REDUCE(orb%p(:,:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%p(:,:,:,:,jspin), 1)
       ENDIF
       CALL MPI_REDUCE(orb%m(:,:,:,:,jspin),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, orb%m(:,:,:,:,jspin), 1)
       ENDIF
       DEALLOCATE (c_b)

    ENDIF

    ! -> Collect the noco staff: 
    IF ( noco%l_noco .AND. jspin.EQ.1 ) THEN

       n = stars%ng3
       ALLOCATE(c_b(n))
       CALL MPI_REDUCE(den%pw(:,3),c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, den%pw(:,3), 1)
       ENDIF
       DEALLOCATE (c_b)
       !
       IF (input%film) THEN

          n = vacuum%nmzxy*(oneD%odi%n2d-1)*2
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(den%vacxy(:,:,:,3),c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, den%vacxy(:,:,:,3), 1)
          ENDIF
          DEALLOCATE (c_b)
          !
          n = vacuum%nmz*2*2
          ALLOCATE(r_b(n))
          CALL MPI_REDUCE(den%vacz(:,:,3:4),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, r_b, 1, den%vacz(:,:,3:4), 1)
          ENDIF
          DEALLOCATE (r_b)

       ENDIF ! input%film


       IF (noco%l_mperp) THEN

          ! -->     for (spin)-off diagonal part of muffin-tin
          n = (atoms%lmaxd+1) * atoms%ntype
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(denCoeffsOffdiag%uu21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%uu21(:,:), 1)
          ENDIF
          CALL MPI_REDUCE(denCoeffsOffdiag%ud21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%ud21(:,:), 1)
          ENDIF
          CALL MPI_REDUCE(denCoeffsOffdiag%du21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%du21(:,:), 1)
          ENDIF
          CALL MPI_REDUCE(denCoeffsOffdiag%dd21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%dd21(:,:), 1)
          ENDIF
          DEALLOCATE (c_b)

          ! -->     lo,u coeff's:
          n = atoms%nlod * atoms%ntype
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(denCoeffsOffdiag%uulo21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%uulo21(:,:), 1)
          ENDIF
          CALL MPI_REDUCE(denCoeffsOffdiag%ulou21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%ulou21(:,:), 1)
          ENDIF
          CALL MPI_REDUCE(denCoeffsOffdiag%dulo21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%dulo21(:,:), 1)
          ENDIF
          CALL MPI_REDUCE(denCoeffsOffdiag%ulod21(:,:),c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%ulod21(:,:), 1)
          ENDIF
          DEALLOCATE (c_b)

          ! -->     lo,lo' coeff's:
          n = atoms%nlod*atoms%nlod*atoms%ntype
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(denCoeffsOffdiag%uloulop21,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%uloulop21, 1)
          ENDIF
          DEALLOCATE (c_b)

          IF (denCoeffsOffdiag%l_fmpl) THEN

             !-->        Full magnetization plots: Collect uunmt21, etc.
             n = (atoms%lmaxd+1)**2 *sphhar%nlhd*atoms%ntype
             ALLOCATE(c_b(n))
             CALL MPI_REDUCE(denCoeffsOffdiag%uunmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%uunmt21, 1)
             ENDIF
             CALL MPI_REDUCE(denCoeffsOffdiag%udnmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%udnmt21, 1)
             ENDIF
             CALL MPI_REDUCE(denCoeffsOffdiag%dunmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%dunmt21, 1)
             ENDIF
             CALL MPI_REDUCE(denCoeffsOffdiag%ddnmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, denCoeffsOffdiag%ddnmt21, 1)
             ENDIF
             DEALLOCATE (c_b)

          ENDIF ! fmpl
       ENDIF  ! mperp
    ENDIF   ! noco

    !+lda+U
    IF ( atoms%n_u.GT.0 ) THEN
       n = 49*atoms%n_u 
       ALLOCATE(c_b(n))
       CALL MPI_REDUCE(n_mmp,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, n_mmp, 1)
       ENDIF
       DEALLOCATE (c_b)
    ENDIF
    !-lda+U

    CALL timestop("mpi_col_den")

  END SUBROUTINE mpi_col_den
END MODULE m_mpi_col_den
