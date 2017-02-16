!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_fjgj
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_fjgj(input,atoms,ispin,cell,lapw,usdus,fj,gj)
!Calculate the fj&gj array which contain the part of the A,B matching coeff. depending on the
!radial functions at the MT boundary as contained in usdus
    USE m_constants, ONLY : fpi_const
    USE m_sphbes
    USE m_dsphbs
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)    :: input
    !TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin 
    !INTEGER, INTENT (IN) :: nintsp
    !LOGICAL,INTENT(IN)   :: l_socfirst

    REAL,INTENT(OUT)     :: fj(:,0:,:,:),gj(:,0:,:,:)
    !     ..
    !     .. Local Scalars ..
    REAL con1,ff,gg,gs,ws

    INTEGER k,l,lo,n,iintsp
    !     .. Local Arrays ..
    REAL gb(0:atoms%lmaxd), fb(0:atoms%lmaxd)
    LOGICAL apw(0:atoms%lmaxd)
    !     ..
    con1 = fpi_const/sqrt(cell%omtil)
    iintsp = ispin
       !$OMP parallel do DEFAULT(SHARED)&
       !$OMP PRIVATE(n,l,apw,lo,k,gs,fb,gb,ws,ff,gg)
       DO n = 1,atoms%ntype
          DO l = 0,atoms%lmax(n)
             apw(l)=any(atoms%l_dulo(:atoms%nlo(n),n))
             IF ((input%l_useapw).AND.(atoms%lapw_l(n).GE.l)) apw(l) = .false.

          ENDDO
          DO lo = 1,atoms%nlo(n)
             IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n)) = .true.
          ENDDO

          DO k = 1,lapw%nv(iintsp)
             gs = lapw%rk(k,iintsp)*atoms%rmt(n)
             CALL sphbes(atoms%lmax(n),gs, fb)
             CALL dsphbs(atoms%lmax(n),gs,fb, gb)
             DO l = 0,atoms%lmax(n)
                !---> set up wronskians for the matching conditions for each ntype
                ws = con1/(usdus%uds(l,n,ispin)*usdus%dus(l,n,ispin)&
                     - usdus%us(l,n,ispin)*usdus%duds(l,n,ispin))
                ff = fb(l)
                gg = lapw%rk(k,iintsp)*gb(l)
                IF ( apw(l) ) THEN
                   fj(k,l,n,iintsp) = 1.0*con1 * ff / usdus%us(l,n,ispin)
                   gj(k,l,n,iintsp) = 0.0d0
                ELSE
                   !IF (noco%l_constr.or.l_socfirst) THEN
                      !--->                 in a constrained calculation fj and gj are needed
                      !--->                 both local spin directions (ispin) at the same time
                   !   fj(k,l,n,ispin) = ws * ( usdus%uds(l,n,ispin)*gg - usdus%duds(l,n,ispin)*ff )
                   !   gj(k,l,n,ispin) = ws * ( usdus%dus(l,n,ispin)*ff - usdus%us(l,n,ispin)*gg )
                   !ELSE
                      !--->                 in a spin-spiral calculation fj and gj are needed
                      !--->                 both interstitial spin directions at the same time
                      !--->                 In all other cases iintsp runs from 1 to 1.
                      fj(k,l,n,iintsp) = ws * ( usdus%uds(l,n,ispin)*gg - usdus%duds(l,n,ispin)*ff )
                      gj(k,l,n,iintsp) = ws * ( usdus%dus(l,n,ispin)*ff - usdus%us(l,n,ispin)*gg )
                   !ENDIF
                ENDIF
             ENDDO
          ENDDO ! k = 1, lapw%nv
       ENDDO    ! n = 1,atoms%ntype
       !$OMP end parallel do

    !ENDDO       ! iintsp = 1,nintsp
  

    RETURN
  END SUBROUTINE hsmt_fjgj
END MODULE m_hsmt_fjgj
