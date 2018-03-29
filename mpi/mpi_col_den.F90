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
  SUBROUTINE mpi_col_den(mpi,sphhar,atoms,oneD,stars,vacuum,&
       input, noco,l_fmpl,jspin,llpd,rhtxy,rht,qpw,ener,&
       sqal,results,svac,pvac,uu,dd,du,uunmt,ddnmt,udnmt,dunmt,sqlo,&
       aclo,bclo,cclo,acnmt,bcnmt,ccnmt,enerlo,orb,mt21,lo21,uloulop21,&
       uunmt21,ddnmt21,udnmt21,dunmt21,den,n_mmp)
    !
#include"cpp_double.h"
    USE m_types
    USE m_constants
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
    INTEGER, INTENT (IN) :: jspin,llpd
    LOGICAL, INTENT (IN) :: l_fmpl
    ! ..
    ! ..  Array Arguments ..
    COMPLEX, INTENT (INOUT) :: qpw(stars%ng3)
    COMPLEX, INTENT (INOUT) :: rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2)
    REAL,    INTENT (INOUT) :: rht(vacuum%nmzd,2) 
    REAL,    INTENT (INOUT) :: ener(0:3,atoms%ntype),sqal(0:3,atoms%ntype)
    REAL,    INTENT (INOUT) :: svac(2),pvac(2)
    REAL,  INTENT (INOUT) :: dd(0:atoms%lmaxd,atoms%ntype)
    REAL,  INTENT (INOUT) :: du(0:atoms%lmaxd,atoms%ntype)
    REAL,  INTENT (INOUT) :: uu(0:atoms%lmaxd,atoms%ntype)
    REAL,  INTENT (INOUT) :: ddnmt(0:llpd,sphhar%nlhd,atoms%ntype)
    REAL,  INTENT (INOUT) :: dunmt(0:llpd,sphhar%nlhd,atoms%ntype)
    REAL,  INTENT (INOUT) :: udnmt(0:llpd,sphhar%nlhd,atoms%ntype)
    REAL,  INTENT (INOUT) :: uunmt(0:llpd,sphhar%nlhd,atoms%ntype)
    REAL,  INTENT (INOUT) :: sqlo(atoms%nlod,atoms%ntype),enerlo(atoms%nlod,atoms%ntype)
    REAL,  INTENT (INOUT) :: aclo(atoms%nlod,atoms%ntype),bclo(atoms%nlod,atoms%ntype)
    REAL,  INTENT (INOUT) :: cclo(atoms%nlod,atoms%nlod,atoms%ntype)
    REAL,  INTENT (INOUT) :: acnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntype)
    REAL,  INTENT (INOUT) :: bcnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntype)
    REAL,  INTENT (INOUT) :: ccnmt(atoms%nlod,atoms%nlod,sphhar%nlhd,atoms%ntype)
    COMPLEX,INTENT(INOUT) :: ddnmt21((atoms%lmaxd+1)**2  )
    COMPLEX,INTENT(INOUT) :: dunmt21((atoms%lmaxd+1)**2  )
    COMPLEX,INTENT(INOUT) :: udnmt21((atoms%lmaxd+1)**2  )
    COMPLEX,INTENT(INOUT) :: uunmt21((atoms%lmaxd+1)**2  )
    COMPLEX,INTENT(INOUT) :: uloulop21(atoms%nlod,atoms%nlod,atoms%ntype)
    COMPLEX,INTENT(INOUT) :: n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u)
    TYPE (t_orb),  INTENT (INOUT) :: orb
    TYPE (t_mt21), INTENT (INOUT) :: mt21(0:atoms%lmaxd,atoms%ntype)
    TYPE (t_lo21), INTENT (INOUT) :: lo21(atoms%nlod,atoms%ntype)
    ! ..
    ! ..  Local Scalars ..
    INTEGER :: n
    ! ..
    ! ..  Local Arrays ..
    INTEGER :: ierr(3)
    COMPLEX, ALLOCATABLE :: c_b(:)
    REAL,    ALLOCATABLE :: r_b(:)
    ! ..
    ! ..  External Subroutines
    EXTERNAL CPP_BLAS_scopy,CPP_BLAS_ccopy,MPI_REDUCE

    !
    ! -> Collect qpw()
    !
    n = stars%ng3
    ALLOCATE(c_b(n))
    CALL MPI_REDUCE(qpw,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_ccopy(n, c_b, 1, qpw, 1)
    ENDIF
    DEALLOCATE (c_b)
    !
    ! -> Collect rhtxy()
    !
    IF (input%film) THEN

       n = vacuum%nmzxyd*(oneD%odi%n2d-1)*2
       ALLOCATE(c_b(n))
       CALL MPI_REDUCE(rhtxy,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_ccopy(n, c_b, 1, rhtxy, 1)
       ENDIF
       DEALLOCATE (c_b)
       !
       ! -> Collect rht()
       !
       n = vacuum%nmzd*2
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(rht,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, rht, 1)
       ENDIF
       DEALLOCATE (r_b)

    ENDIF
    !
    ! -> Collect uu(),ud() and dd()
    !
    n = (atoms%lmaxd+1)*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(uu,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, uu, 1)
    ENDIF
    CALL MPI_REDUCE(du,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, du, 1)
    ENDIF
    CALL MPI_REDUCE(dd,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, dd, 1)
    ENDIF
    DEALLOCATE (r_b)
    !
    !--> Collect uunmt,udnmt,dunmt,ddnmt
    !
    n = (llpd+1)*sphhar%nlhd*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(uunmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, uunmt, 1)
    ENDIF
    CALL MPI_REDUCE(udnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, udnmt, 1)
    ENDIF
    CALL MPI_REDUCE(dunmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, dunmt, 1)
    ENDIF
    CALL MPI_REDUCE(ddnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, ddnmt, 1)
    ENDIF
    DEALLOCATE (r_b)
    !
    !--> ener & sqal
    !
    n=4*atoms%ntype
    ALLOCATE(r_b(n))
    CALL MPI_REDUCE(ener,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, ener, 1)
    ENDIF
    CALL MPI_REDUCE(sqal,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
    IF (mpi%irank.EQ.0) THEN
       CALL CPP_BLAS_scopy(n, r_b, 1, sqal, 1)
    ENDIF
    DEALLOCATE (r_b)
    !
    !--> svac & pvac
    !
    IF ( input%film ) THEN

       n=2
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(svac,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, svac, 1)
       ENDIF
       CALL MPI_REDUCE(pvac,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, pvac, 1)
       ENDIF
       DEALLOCATE (r_b)

    ENDIF
    !
    ! -> Collect force
    !   
    IF (input%l_f) THEN

       n=3*atoms%ntype
       ALLOCATE(r_b(n))
       CALL MPI_REDUCE(results%force(1,1,jspin),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, results%force(1,1,jspin), 1)
       ENDIF
       DEALLOCATE (r_b)

    ENDIF
    !
    ! -> Optional the LO-coefficients: aclo,bclo,enerlo,cclo,acnmt,bcnmt,ccnmt
    !
    IF (atoms%nlod.GE.1) THEN

       n=atoms%nlod*atoms%ntype 
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(aclo,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, aclo, 1)
       ENDIF
       CALL MPI_REDUCE(bclo,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, bclo, 1)
       ENDIF
       CALL MPI_REDUCE(enerlo,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, enerlo, 1)
       ENDIF
       CALL MPI_REDUCE(sqlo,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, sqlo, 1)
       ENDIF
       DEALLOCATE (r_b)

       n = atoms%nlod * atoms%nlod * atoms%ntype
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(cclo,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, cclo, 1)
       ENDIF
       DEALLOCATE (r_b)

       n = (atoms%lmaxd+1) * atoms%ntype * atoms%nlod * sphhar%nlhd
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(acnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, acnmt, 1)
       ENDIF
       CALL MPI_REDUCE(bcnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, bcnmt, 1)
       ENDIF
       DEALLOCATE (r_b)

       n = atoms%ntype * sphhar%nlhd * atoms%nlod**2
       ALLOCATE (r_b(n))
       CALL MPI_REDUCE(ccnmt,r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
       IF (mpi%irank.EQ.0) THEN
          CALL CPP_BLAS_scopy(n, r_b, 1, ccnmt, 1)
       ENDIF
       DEALLOCATE (r_b)

    ENDIF
    !
    ! ->  Now the SOC - stuff: orb, orblo and orblo
    !
    IF (noco%l_soc) THEN
       !
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

    !
    ! -> Collect the noco staff: 
    !
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

          n = vacuum%nmzxyd*(oneD%odi%n2d-1)*2
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(den%vacxy(:,:,:,3),c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, den%vacxy(:,:,:,3), 1)
          ENDIF
          DEALLOCATE (c_b)
          !
          n = vacuum%nmzd*2*2
          ALLOCATE(r_b(n))
          CALL MPI_REDUCE(den%vacz(:,:,3:4),r_b,n,CPP_MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, r_b, 1, den%vacz(:,:,3:4), 1)
          ENDIF
          DEALLOCATE (r_b)

       ENDIF ! input%film


       IF (noco%l_mperp) THEN
          !
          ! -->     for (spin)-off diagonal part of muffin-tin
          !
          n = (atoms%lmaxd+1) * atoms%ntype
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(mt21(:,:)%uu,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%uu, 1)
          ENDIF
          CALL MPI_REDUCE(mt21(:,:)%ud,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%ud, 1)
          ENDIF
          CALL MPI_REDUCE(mt21(:,:)%du,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%du, 1)
          ENDIF
          CALL MPI_REDUCE(mt21(:,:)%dd,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, mt21(:,:)%dd, 1)
          ENDIF
          DEALLOCATE (c_b)
          !
          ! -->     lo,u coeff's:
          !
          n = atoms%nlod * atoms%ntype
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(lo21(:,:)%uulo,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%uulo, 1)
          ENDIF
          CALL MPI_REDUCE(lo21(:,:)%ulou,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%ulou, 1)
          ENDIF
          CALL MPI_REDUCE(lo21(:,:)%dulo,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%dulo, 1)
          ENDIF
          CALL MPI_REDUCE(lo21(:,:)%ulod,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, lo21(:,:)%ulod, 1)
          ENDIF
          DEALLOCATE (c_b)
          !
          ! -->     lo,lo' coeff's:
          !
          n = atoms%nlod*atoms%nlod*atoms%ntype
          ALLOCATE(c_b(n))
          CALL MPI_REDUCE(uloulop21,c_b,n,CPP_MPI_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD,ierr)
          IF (mpi%irank.EQ.0) THEN
             CALL CPP_BLAS_ccopy(n, c_b, 1, uloulop21, 1)
          ENDIF
          DEALLOCATE (c_b)

          IF (l_fmpl) THEN
             !
             !-->        Full magnetization plots: Collect uunmt21, etc.
             !
             n = (atoms%lmaxd+1)**2 *sphhar%nlhd*atoms%ntype
             ALLOCATE(c_b(n))
             CALL MPI_REDUCE(uunmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, uunmt21, 1)
             ENDIF
             CALL MPI_REDUCE(udnmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, udnmt21, 1)
             ENDIF
             CALL MPI_REDUCE(dunmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, dunmt21, 1)
             ENDIF
             CALL MPI_REDUCE(ddnmt21,c_b,n,CPP_MPI_COMPLEX,MPI_SUM,0 ,MPI_COMM_WORLD,ierr)
             IF (mpi%irank.EQ.0) THEN
                CALL CPP_BLAS_ccopy(n, c_b, 1, ddnmt21, 1)
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

    RETURN
  END SUBROUTINE mpi_col_den
END MODULE m_mpi_col_den
