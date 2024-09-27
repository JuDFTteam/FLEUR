!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_fjgj
  USE m_juDFT
  IMPLICIT NONE

  PRIVATE
  TYPE t_fjgj
    REAL,ALLOCATABLE    :: fj(:,:,:,:),gj(:,:,:,:)
  CONTAINS
    procedure :: alloc
    procedure :: calculate => hsmt_fjgj_cpu
  END TYPE
  PUBLIC t_fjgj

CONTAINS
  subroutine alloc(fjgj,nvd,lmaxd,isp,noco)
    USE m_types
    CLASS(t_fjgj),INTENT(OUT) :: fjgj
    INTEGER,INTENT(IN)        :: nvd,lmaxd,isp
    TYPE(t_noco),INTENT(IN)   :: noco

    ALLOCATE(fjgj%fj(nvd,0:lmaxd,merge(1,isp,noco%l_noco):merge(2,isp,noco%l_noco),MERGE(2,1,noco%l_ss)))
    ALLOCATE(fjgj%gj(nvd,0:lmaxd,merge(1,isp,noco%l_noco):merge(2,isp,noco%l_noco),MERGE(2,1,noco%l_ss)))

    fjgj%fj = 0.0
    fjgj%gj = 0.0

  end subroutine

  SUBROUTINE hsmt_fjgj_cpu(fjgj,input,atoms,cell,lapw,noco,usdus,n,ispin)
    !Calculate the fj&gj array which contain the part of the A,B matching coeff. depending on the
    !radial functions at the MT boundary as contained in usdus
    USE m_constants, ONLY : fpi_const
    USE m_sphbes
    USE m_dsphbs
    USE m_types
    IMPLICIT NONE
    CLASS(t_fjgj),INTENT(INOUT) :: fjgj
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin,n

    !     ..
    !     .. Local Scalars ..
    REAL con1,ff,gg,gs,ws

    INTEGER k,l,lo,intspin,jspin, jspinStart, jSpinEnd
    LOGICAL l_socfirst
    !     .. Local Arrays ..
    REAL gb(0:atoms%lmaxd), fb(0:atoms%lmaxd)
    LOGICAL apw(0:atoms%lmaxd)
    !     ..
    l_socfirst = noco%l_soc .AND. noco%l_noco .AND. (.NOT. noco%l_ss)
    con1 = fpi_const/SQRT(cell%omtil)
    DO l = 0,atoms%lmax(n)
       apw(l)=ANY(atoms%l_dulo(:atoms%nlo(n),n))
       IF ((input%l_useapw).AND.(atoms%lapw_l(n).GE.l)) apw(l) = .FALSE.
    ENDDO
    DO lo = 1,atoms%nlo(n)
       IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n)) = .TRUE.
    ENDDO

    jspinStart = ispin
    jspinEnd = ispin
    IF (any(noco%l_constrained).or.l_socfirst.OR.any(noco%l_unrestrictMT).OR.any(noco%l_spinoffd_ldau)) THEN
       jspinStart = 1
       jspinEnd = input%jspins
    END IF

    DO intspin=1,MERGE(2,1,noco%l_ss)
#ifndef _OPENACC
       !$OMP PARALLEL DO DEFAULT(NONE) &
       !$OMP PRIVATE(l,gs,fb,gb,ws,ff,gg,jspin)&
       !$OMP SHARED(lapw,atoms,con1,usdus,l_socfirst,noco,input)&
       !$OMP SHARED(fjgj,intspin,n,ispin,apw,jspinStart,jspinEnd)
#else
       !$acc kernels present(fjgj,fjgj%fj,fjgj%gj) &
       !$acc &copyin(lapw,lapw%rk,lapw%nv,atoms,atoms%rmt,atoms%lmax,apw,con1)&
       !$acc &copyin(usdus,usdus%dus,usdus%uds,usdus%us,usdus%duds)&
       !$acc &copyin(intspin,n)
       !$acc loop gang private(k,gs,fb,gb)!,ws,ff,gg,jspin)
#endif
       DO k = 1,lapw%nv(intspin)
          gs = lapw%rk(k,intspin)*atoms%rmt(n)
          CALL sphbes(atoms%lmax(n),gs, fb)
          CALL dsphbs(atoms%lmax(n),gs,fb, gb)
!          !$OMP SIMD PRIVATE(ws,ff,gg)
          !$acc loop vector PRIVATE(l,ws,ff,gg,jspin)
          DO l = 0,atoms%lmax(n)
             ff = fb(l)
             gg = lapw%rk(k,intspin)*gb(l)
             DO jspin = jspinStart, jspinEnd
                IF ( apw(l) ) THEN
                   fjgj%fj(k,l,jspin,intspin) = 1.0*con1 * ff / usdus%us(l,n,jspin)
                   fjgj%gj(k,l,jspin,intspin) = 0.0
                ELSE
                   ws = con1/(usdus%uds(l,n,jspin)*usdus%dus(l,n,jspin)- usdus%us(l,n,jspin)*usdus%duds(l,n,jspin))
                   fjgj%fj(k,l,jspin,intspin) = ws * ( usdus%uds(l,n,jspin)*gg - usdus%duds(l,n,jspin)*ff )
                   fjgj%gj(k,l,jspin,intspin) = ws * ( usdus%dus(l,n,jspin)*ff - usdus%us(l,n,jspin)*gg )
                ENDIF
             END DO
          ENDDO
          !$acc end loop
!          !$OMP END SIMD
       ENDDO ! k = 1, lapw%nv
#ifdef _OPENACC
       !$acc end loop
       !$acc end kernels
#else
       !$OMP END PARALLEL DO
#endif
    ENDDO
    RETURN
  END SUBROUTINE hsmt_fjgj_cpu
END MODULE m_hsmt_fjgj
