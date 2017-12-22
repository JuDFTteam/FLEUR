!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_fjgj
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_fjgj(input,atoms,cell,lapw,noco,usdus,n,ispin,fj,gj)
    !Calculate the fj&gj array which contain the part of the A,B matching coeff. depending on the
    !radial functions at the MT boundary as contained in usdus
    USE m_constants, ONLY : fpi_const
    USE m_sphbes
    USE m_dsphbs
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin,n
  
    REAL,INTENT(OUT)     :: fj(:,0:,:),gj(:,0:,:)
    !     ..
    !     .. Local Scalars ..
    REAL con1,ff,gg,gs,ws

    INTEGER k,l,lo,intspin
    !     .. Local Arrays ..
    REAL gb(0:atoms%lmaxd), fb(0:atoms%lmaxd)
    LOGICAL apw(0:atoms%lmaxd)
    !     ..
    con1 = fpi_const/SQRT(cell%omtil)
    DO l = 0,atoms%lmax(n)
       apw(l)=ANY(atoms%l_dulo(:atoms%nlo(n),n))
       IF ((input%l_useapw).AND.(atoms%lapw_l(n).GE.l)) apw(l) = .FALSE.
    ENDDO
    DO lo = 1,atoms%nlo(n)
       IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n)) = .TRUE.
    ENDDO
    DO intspin=1,MERGE(2,1,noco%l_ss)
       !$OMP PARALLEL DO DEFAULT(NONE) &
       !$OMP PRIVATE(l,gs,fb,gb,ws,ff,gg)&
       !$OMP SHARED(lapw,atoms,con1,usdus)&
       !$OMP SHARED(fj,gj,intspin,n,ispin,apw)
       DO k = 1,lapw%nv(intspin)
          gs = lapw%rk(k,intspin)*atoms%rmt(n)
          CALL sphbes(atoms%lmax(n),gs, fb)
          CALL dsphbs(atoms%lmax(n),gs,fb, gb)
          !$OMP SIMD PRIVATE(ws,ff,gg)
          DO l = 0,atoms%lmax(n)
             !---> set up wronskians for the matching conditions for each ntype
             ws = con1/(usdus%uds(l,n,ispin)*usdus%dus(l,n,ispin)&
                  - usdus%us(l,n,ispin)*usdus%duds(l,n,ispin))
             ff = fb(l)
             gg = lapw%rk(k,intspin)*gb(l)
             IF ( apw(l) ) THEN
                fj(k,l,intspin) = 1.0*con1 * ff / usdus%us(l,n,ispin)
                gj(k,l,intspin) = 0.0d0
             ELSE
                fj(k,l,intspin) = ws * ( usdus%uds(l,n,ispin)*gg - usdus%duds(l,n,ispin)*ff )
                gj(k,l,intspin) = ws * ( usdus%dus(l,n,ispin)*ff - usdus%us(l,n,ispin)*gg )
                !ENDIF
             ENDIF
          ENDDO
          !$OMP END SIMD
       ENDDO ! k = 1, lapw%nv
       !$OMP END PARALLEL DO
    ENDDO
    RETURN
  END SUBROUTINE hsmt_fjgj
END MODULE m_hsmt_fjgj
