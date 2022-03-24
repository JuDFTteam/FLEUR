!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsvac
   USE m_juDFT
CONTAINS
   !-----------------------------------------------------------------------------
   ! Calculate the vacuum contribution to the Hamiltonian and Overlap matrix
   !-----------------------------------------------------------------------------
   SUBROUTINE hsvac(vacuum, stars, fmpi, jsp, input, v, evac, cell, &
                  & lapw, sym, noco, nococonv, hmat, smat)

      USE m_vacfun
      USE m_types

      IMPLICIT NONE

      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_vacuum),INTENT(IN)     :: vacuum
      TYPE(t_noco),INTENT(IN)       :: noco
      TYPE(t_nococonv),INTENT(IN)       :: nococonv
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_stars),INTENT(IN)      :: stars
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_lapw),INTENT(IN)       :: lapw
      TYPE(t_mpi),INTENT(IN)        :: fmpi
      TYPE(t_potden),INTENT(IN)     :: v
      CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:),smat(:,:)
      !     ..
      !     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jsp
      !     ..
      !     .. Array Arguments ..
      REAL,    INTENT (IN) :: evac(2,input%jspins)
      !     ..
      !     .. Local Scalars ..
      COMPLEX :: hij,sij,apw_lo,c_1
      REAL    :: d2,gz,sign,th,wronk
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
      REAL    :: uz(lapw%dim_nv2d(),input%jspins)

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
                         & sym, cell, ivac, evac, lapw%bkpt, v%vacxy, v%vacz, kvac, nv2, &
                         & tuuv, tddv, tudv, tduv, uz, duz, udz, dudz, ddnv, wronk)
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
               END DO

               !---> update hamiltonian and overlap matrices
               IF (iSpin==iSpinPr) THEN
                  DO ikG = fmpi%n_rank+1,lapw%nv(iSpinPr),fmpi%n_size
                     ikG0=(ikG-1)/fmpi%n_size+1 !local column index
                     ikG2 = map2(ikG,iSpinPr)
                     DO ikGPr = 1,ikG - 1 !TODO check noco case
                        !---> overlap: only  (g-g') parallel=0       '
                        IF (map2(ikGPr,iSpin).EQ.ikG2) THEN
                           sij = CONJG(a(ikG,iSpinPr))*a(ikGPr,iSpinPr) + &
                                 CONJG(b(ikG,iSpinPr))*b(ikGPr,iSpinPr)*ddnv(ikG2,iSpinPr)
                           !+APW_LO
                           IF (input%l_useapw) THEN
                              apw_lo = CONJG(a(ikG,iSpin)*  uz(ikG2,iSpin) + b(ikG,iSpin)* udz(ikG2,iSpin) ) &
                                          * (a(ikGPr,iSpin)* duz(ikG2,iSpin) + b(ikGPr,iSpin)*dudz(ikG2,iSpin) )&
                                     +      (a(ikGPr,iSpin)*  uz(ikG2,iSpin) + b(ikGPr,iSpin)* udz(ikG2,iSpin) ) &
                                     * CONJG(a(ikG,iSpin)* duz(ikG2,iSpin) + b(ikG,iSpin)*dudz(ikG2,iSpin) )
                              ! IF (i.lt.10) write (3,'(2i4,2f20.10)') i,j,apw_lo
                              IF (hmat(1,1)%l_real) THEN
                                 hmat(igSpin,igSpin)%data_r(ikGPr,ikG0) = hmat(igSpin,igSpin)%data_r(ikGPr,ikG0) + 0.25 * REAL(apw_lo)
                              ELSE
                                 hmat(igSpin,igSpin)%data_c(ikGPr,ikG0) = hmat(igSpin,igSpin)%data_c(ikGPr,ikG0) + 0.25 * apw_lo
                              END IF
                           END IF

                           !Overlap matrix
                           IF (hmat(1,1)%l_real) THEN
                              smat(igSpin,igSpin)%data_r(ikGPr,ikG0) = smat(igSpin,igSpin)%data_r(ikGPr,ikG0) + REAL(sij)
                           ELSE
                              smat(igSpin,igSpin)%data_c(ikGPr,ikG0) = smat(igSpin,igSpin)%data_c(ikGPr,ikG0) + sij
                           END IF
                        END IF
                     END DO

                     !Diagonal term of Overlapp matrix, Hamiltonian later
                     sij = CONJG(a(ikG,iSpinPr))*a(ikG,iSpinPr) + CONJG(b(ikG,iSpinPr))*b(ikG,iSpinPr)*ddnv(ikG2,iSpinPr)
                     IF (hmat(1,1)%l_real) THEN
                        smat(igSpin,igSpin)%data_r(ikGPr,ikG0) = smat(igSpin,igSpin)%data_r(ikGPr,ikG0) + REAL(sij)
                     ELSE
                        smat(igSpin,igSpin)%data_c(ikGPr,ikG0) = smat(igSpin,igSpin)%data_c(ikGPr,ikG0) + sij
                     END IF
                  END DO
               END IF

               !--->    hamiltonian update
               DO  ikG = fmpi%n_rank+1,lapw%nv(iSpin),fmpi%n_size
                  ikG0=(ikG-1)/fmpi%n_size+1 !local column index
                  ikG2 = map2(ikG,iSpin)
                  DO ikGPr = 1,MERGE(ikG,lapw%nv(iSpinPr),iSpin==iSpinPr)
                     ikGPr2 = map2(ikGPr,iSpinPr)
                     hij = CONJG(a(ikG,iSpin))* (tuuv(ikG2,ikGPr2)*a(ikGPr,iSpinPr) +tudv(ikG2,ikGPr2)*b(ikGPr,iSpinPr))&
                         + CONJG(b(ikG,iSpin))* (tddv(ikG2,ikGPr2)*b(ikGPr,iSpinPr) +tduv(ikG2,ikGPr2)*a(ikGPr,iSpinPr))
                     IF (hmat(1,1)%l_real) THEN
                        hmat(igSpinPr,igSpin)%data_r(ikGPr,ikG0) = hmat(igSpinPr,igSpin)%data_r(ikGPr,ikG0) + REAL(hij)
                     ELSE
                        hmat(igSpinPr,igSpin)%data_c(ikGPr,ikG0) = hmat(igSpinPr,igSpin)%data_c(ikGPr,ikG0) + hij
                     END IF
                  END DO
               END DO
               !--->    end of loop over different parts of the potential matrix
            END DO
            !---> end of loop over vacua
         END DO
      END DO
   END SUBROUTINE hsvac
END MODULE m_hsvac
