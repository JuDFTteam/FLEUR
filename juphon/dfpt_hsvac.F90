!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_hsvac
   USE m_juDFT
CONTAINS
   !-----------------------------------------------------------------------------
   ! Calculate the vacuum contribution to the Hamiltonian and Overlap matrix
   !-----------------------------------------------------------------------------
   SUBROUTINE dfpt_hsvac(vacuum, stars, fmpi, jsp, input, v, v1, evac, cell, &
                  & lapwq, lapw, noco, nococonv, hmat)

      USE m_vacfun
      USE m_types

      IMPLICIT NONE

      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_vacuum),INTENT(IN)     :: vacuum
      TYPE(t_noco),INTENT(IN)       :: noco
      TYPE(t_nococonv),INTENT(IN)       :: nococonv
      TYPE(t_stars),INTENT(IN)      :: stars
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_lapw),INTENT(IN)       :: lapw, lapwq
      TYPE(t_mpi),INTENT(IN)        :: fmpi
      TYPE(t_potden),INTENT(IN)     :: v, v1
      CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)!,smat(:,:)
      !     ..
      !     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jsp
      !     ..
      !     .. Array Arguments ..
      REAL,    INTENT (IN) :: evac(2,input%jspins)
      !     ..
      !     .. Local Scalars ..
      COMPLEX :: hij,sij,apw_lo,c_1
      REAL    :: d2,gz,sign,th,wronk,wronkq
      INTEGER :: ikG,ikG2,ikGPr,ikGPr2,jspin,ikG0
      INTEGER :: ivac,igSpin,igSpinPr
      INTEGER :: iSpin,iSpinPr
      INTEGER :: nc
      !     ..
      !     .. Local Arrays ..
      INTEGER :: nv2(input%jspins)
      INTEGER :: kvac(2,lapw%dim_nv2d(),input%jspins)
      INTEGER :: map2(lapw%dim_nvd(),input%jspins)
      COMPLEX :: tddv(lapw%dim_nv2d(),lapw%dim_nv2d()),tduv(lapw%dim_nv2d(),lapw%dim_nv2d())
      COMPLEX :: tudv(lapw%dim_nv2d(),lapw%dim_nv2d()),tuuv(lapw%dim_nv2d(),lapw%dim_nv2d())
      COMPLEX :: a(lapw%dim_nvd(),input%jspins),b(lapw%dim_nvd(),input%jspins)
      REAL    :: ddnv(lapw%dim_nv2d(),input%jspins),dudz(lapw%dim_nv2d(),input%jspins)
      REAL    :: duz(lapw%dim_nv2d(),input%jspins), udz(lapw%dim_nv2d(),input%jspins)
      REAL    :: uz(lapw%dim_nv2d(),input%jspins), dudzq(lapw%dim_nv2d(),input%jspins)
      REAL    :: duzq(lapw%dim_nv2d(),input%jspins), udzq(lapw%dim_nv2d(),input%jspins)
      REAL    :: uzq(lapw%dim_nv2d(),input%jspins)
      COMPLEX :: aq(lapw%dim_nvd(),input%jspins),bq(lapw%dim_nvd(),input%jspins)
      INTEGER :: map2q(lapw%dim_nvd(),input%jspins)
      INTEGER :: nv2q(input%jspins)
      INTEGER :: kvacq(2,lapw%dim_nv2d(),input%jspins)

      d2 = SQRT(cell%omtil/cell%area)

      !---> set up mapping function from 3d-->2d lapws
      DO jspin = 1,input%jspins
         nv2(jspin) = 0
         k_loop:DO ikG = 1, lapw%nv(jspin)
            DO ikG2 = 1, nv2(jspin)
               IF (all(lapw%gvec(1:2,ikG,jspin)==kvac(1:2,ikG2,jspin))) THEN
                  map2(ikG,jspin) = ikG2
                  CYCLE k_loop
               END IF
            END DO
            nv2(jspin) = nv2(jspin) + 1
            IF (nv2(jspin)>lapw%dim_nv2d())  CALL juDFT_error("hsvac:lapw%dim_nv2d()",calledby ="hsvac")
            kvac(1:2,nv2(jspin),jspin) = lapw%gvec(1:2,ikG,jspin)
            map2(ikG,jspin) = nv2(jspin)
         END DO k_loop
      END DO

      DO jspin = 1,input%jspins
         nv2q(jspin) = 0
         k_loopq:DO ikG = 1, lapwq%nv(jspin)
            DO ikG2 = 1, nv2q(jspin)
               IF (all(lapwq%gvec(1:2,ikG,jspin)==kvacq(1:2,ikG2,jspin))) THEN
                  map2q(ikG,jspin) = ikG2
                  CYCLE k_loopq
               END IF
            END DO
            nv2q(jspin) = nv2q(jspin) + 1
            IF (nv2q(jspin)>lapw%dim_nv2d()) CALL juDFT_error("hsvac:lapw%dim_nv2d()",calledby ="hsvac")
            kvacq(1:2,nv2q(jspin),jspin) = lapwq%gvec(1:2,ikG,jspin)
            map2q(ikG,jspin) = nv2q(jspin)
         END DO k_loopq
      END DO

      !---> loop over the two vacuua (1: upper; 2: lower)
      DO ivac = 1,2
         sign = 3. - 2.*ivac !+/- 1
         DO iSpin=MERGE(1,jsp,noco%l_noco),MERGE(2,jsp,noco%l_noco) !loop over global spin
            igSpin=MIN(SIZE(hmat,1),iSpin) !in colinear case igSpin=1
            DO iSpinPr=MERGE(1,jsp,noco%l_noco),MERGE(2,jsp,noco%l_noco) !loop over global spin
               igSpinPr=MIN(SIZE(hmat,1),iSpinPr) !in colinear case igSpinPr=1
               !---> get the wavefunctions and set up the tuuv, etc matrices
               CALL timestart("vacfun")
               CALL vacfun(fmpi, vacuum, stars, input, nococonv, iSpin, iSpinPr, &
                         & cell, ivac, evac, lapw%bkpt, v%vac(:vacuum%nmzxyd,2:,:,:), v%vac(:,1,:,:), kvac, nv2, &
                         & tuuv, tddv, tudv, tduv, uz, duz, udz, dudz, ddnv, wronk,&
                         & lapwq%bkpt+lapwq%qphon, v1%vac(:vacuum%nmzxyd,2:,:,:), v1%vac(:,1,:,:), kvacq, nv2q, uzq, duzq, udzq, dudzq, wronkq)
               CALL timestop("vacfun")

               !---> generate a and b coeffficients
               DO jspin = MIN(iSpin,iSpinPr),MAX(iSpin,iSpinPr)
                  DO ikG = 1,lapw%nv(jspin)
                     gz = sign*cell%bmat(3,3)*lapw%k3(ikG,jspin)
                     ikG2 = map2(ikG,jspin)
                     th = gz*cell%z1
                     c_1 = CMPLX( COS(th), SIN(th) )/ (d2*wronk)
                     a(ikG,jspin) = - c_1 * CMPLX(dudz(ikG2,jspin), gz*udz(ikG2,jspin) )
                     b(ikG,jspin) =   c_1 * CMPLX(duz(ikG2,jspin), gz* uz(ikG2,jspin) )
                  END DO

                  DO ikG = 1,lapwq%nv(jspin)
                     gz = sign*cell%bmat(3,3)*lapwq%k3(ikG,jspin)
                     ikG2 = map2q(ikG,jspin)
                     th = gz*cell%z1
                     c_1 = CMPLX( COS(th), SIN(th) )/ (d2*wronkq)
                     aq(ikG,jspin) = - c_1 * CMPLX(dudzq(ikG2,jspin), gz*udzq(ikG2,jspin) )
                     bq(ikG,jspin) =   c_1 * CMPLX(duzq(ikG2,jspin), gz* uzq(ikG2,jspin) )
                  END DO
               END DO

               !---> update hamiltonian and overlap matrices
               !!IF (iSpinPr==iSpin) THEN
               !!   DO ikG = fmpi%n_rank + 1, lapw%nv(iSpin), fmpi%n_size
               !!      ikG0 = (ikG-1)/fmpi%n_size + 1 !local column index
               !!      ikG2 = map2(ikG,iSpin)
               !!      DO ikGPr = 1, ikG - 1 !TODO check noco case
               !!         !---> overlap: only  (g-g') parallel=0       '
               !!         IF (map2(ikGPr, iSpin).EQ.ikG2) THEN
               !!            sij = CONJG(a(ikGPr,iSpin))*a(ikG,iSpin) + &
               !!                  CONJG(b(ikGPr,iSpin))*b(ikG,iSpin)*ddnv(ikG2,iSpin)
               !!            !+APW_LO
               !!            IF (input%l_useapw) THEN
               !!               apw_lo =      (a(ikG,iSpin)   *  uz(ikG2,iSpin) + b(ikG,iSpin)   *  udz(ikG2,iSpin)) &
               !!                     * CONJG(a(ikGPr,iSpin) * duz(ikG2,iSpin) + b(ikGPr,iSpin) * dudz(ikG2,iSpin)) &
               !!                      + CONJG(a(ikGPr,iSpin) *  uz(ikG2,iSpin) + b(ikGPr,iSpin) *  udz(ikG2,iSpin)) &
               !!                      *      (a(ikG,iSpin)   * duz(ikG2,iSpin) + b(ikG,iSpin)   * dudz(ikG2,iSpin))
               !!               ! IF (i.lt.10) write (3,'(2i4,2f20.10)') i,j,apw_lo
               !!               IF (hmat(1,1)%l_real) THEN
               !!                  hmat(igSpin,igSpin)%data_r(ikGPr,ikG0) = hmat(igSpin,igSpin)%data_r(ikGPr,ikG0) + 0.25 * REAL(apw_lo)
               !!               ELSE
               !!                  hmat(igSpin,igSpin)%data_c(ikGPr,ikG0) = hmat(igSpin,igSpin)%data_c(ikGPr,ikG0) + 0.25 * apw_lo
               !!               END IF
               !!            END IF

               !!            !Overlap matrix
               !!            IF (hmat(1,1)%l_real) THEN
               !!               smat(igSpin,igSpin)%data_r(ikGPr,ikG0) = smat(igSpin,igSpin)%data_r(ikGPr,ikG0) + REAL(sij)
               !!            ELSE
               !!               smat(igSpin,igSpin)%data_c(ikGPr,ikG0) = smat(igSpin,igSpin)%data_c(ikGPr,ikG0) + sij
               !!            END IF
               !!        END IF
               !!      END DO

               !!      !Diagonal term of Overlap matrix, Hamiltonian later
               !!      sij = CONJG(a(ikG,iSpin))*a(ikG,iSpin) + CONJG(b(ikG,iSpin))*b(ikG,iSpin)*ddnv(ikG2,iSpin)
               !!      IF (hmat(1,1)%l_real) THEN
               !!         smat(igSpin,igSpin)%data_r(ikGPr,ikG0) = smat(igSpin,igSpin)%data_r(ikGPr,ikG0) + REAL(sij)
               !!      ELSE
               !!         smat(igSpin,igSpin)%data_c(ikGPr,ikG0) = smat(igSpin,igSpin)%data_c(ikGPr,ikG0) + sij
               !!      END IF
               !!   END DO
               !!END IF

               !--->    hamiltonian update
               DO ikG = fmpi%n_rank+1,lapw%nv(iSpin),fmpi%n_size
                  ikG0 = (ikG-1)/fmpi%n_size + 1 !local column index
                  ikG2 = map2(ikG,iSpin)
                  DO ikGPr = 1, lapwq%nv(iSpinPr)
                     ikGPr2 = map2q(ikGPr, iSpinPr)
                     hij = CONJG(aq(ikGPr, iSpinPr)) * tuuv(ikGPr2, ikG2) * a(ikG,iSpin) &
                         + CONJG(bq(ikGPr, iSpinPr)) * tddv(ikGPr2, ikG2) * b(ikG,iSpin) &
                         + CONJG(aq(ikGPr, iSpinPr)) * tudv(ikGPr2, ikG2) * b(ikG,iSpin) &
                         + CONJG(bq(ikGPr, iSpinPr)) * tduv(ikGPr2, ikG2) * a(ikG,iSpin)
                     hmat(igSpinPr,igSpin)%data_c(ikGPr,ikG0) = hmat(igSpinPr,igSpin)%data_c(ikGPr,ikG0) + hij
                  END DO
               END DO
               !--->    end of loop over different parts of the potential matrix
            END DO
            !---> end of loop over vacua
         END DO
      END DO
   END SUBROUTINE dfpt_hsvac
END MODULE m_dfpt_hsvac
