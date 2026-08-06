!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_hsvac
   USE m_juDFT
   implicit none
CONTAINS
   !-----------------------------------------------------------------------------
   ! Calculate the vacuum contribution to the Hamiltonian perturbation
   !-----------------------------------------------------------------------------
   SUBROUTINE dfpt_hsvac(vacuum, stars, fmpi, jsp, input, v, v1, evac, cell, &
                  & lapwq, lapw, noco, nococonv, hmat)

      USE m_vacfun
      USE m_vac_map2
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
      CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)

      INTEGER, INTENT (IN) :: jsp
      REAL,    INTENT (IN) :: evac(2,input%jspins)

      COMPLEX :: hij,c_1
      REAL    :: d2,gz,sign,th,wronk
      INTEGER :: ikG,ikG2,ikGPr,ikGPr2,jspin,ikG0
      INTEGER :: ivac,igSpin,igSpinPr
      INTEGER :: iSpin,iSpinPr

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

      !---> set up mapping function from 3d-->2d lapws, at k and at k+q
      CALL vac_map2(lapw,  input%jspins, nv2,  kvac,  map2)
      CALL vac_map2(lapwq, input%jspins, nv2q, kvacq, map2q)

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
                         & lapwq%bkpt+lapwq%qphon, v1%vac, kvacq, nv2q, uzq, duzq, udzq, dudzq)
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
                     c_1 = CMPLX( COS(th), SIN(th) )/ (d2*wronk)
                     aq(ikG,jspin) = - c_1 * CMPLX(dudzq(ikG2,jspin), gz*udzq(ikG2,jspin) )
                     bq(ikG,jspin) =   c_1 * CMPLX(duzq(ikG2,jspin), gz* uzq(ikG2,jspin) )
                  END DO
               END DO

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
