!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_fjgj
  USE m_juDFT
  IMPLICIT NONE

  INTERFACE hsmt_fjgj
    module procedure hsmt_fjgj_cpu
#ifdef CPP_GPU
    module procedure hsmt_fjgj_gpu
#endif
  END INTERFACE

CONTAINS
#ifdef CPP_GPU

  SUBROUTINE synth_fjgj(nv,ispin,jspins,lmax,lmaxd,apw,l_flag,rk,rmt,con1,uds,dus,us,duds,fj,gj)
  USE m_sphbes
  USE m_dsphbs
  INTEGER, INTENT(IN) :: nv, ispin, jspins, lmax, lmaxd
  LOGICAL, INTENT(IN) :: apw(0:lmaxd), l_flag
  REAL, INTENT(IN) :: rk(:),rmt,con1
  REAL, INTENT(IN) :: uds(0:lmaxd,jspins),dus(0:lmaxd,jspins),us(0:lmaxd,jspins),duds(0:lmaxd,jspins) 
  REAL,INTENT(OUT),MANAGED     :: fj(:,0:,:),gj(:,0:,:)

  REAL gb(0:lmaxd), fb(0:lmaxd)
  REAL ws(jspins)
  REAL ff,gg,gs
  INTEGER k,l,jspin

  DO k = 1,nv
          gs = rk(k)*rmt
          CALL sphbes(lmax,gs,fb)
          CALL dsphbs(lmax,gs,fb,gb)
          DO l = 0,lmax
             !---> set up wronskians for the matching conditions for each ntype
             DO jspin = 1, jspins
                ws(jspin) = con1/(uds(l,jspin)*dus(l,jspin) - us(l,jspin)*duds(l,jspin))
             END DO
             ff = fb(l)
             gg = rk(k)*gb(l)
             IF ( apw(l) ) THEN
                fj(k,l,ispin) = 1.0*con1 * ff / us(l,ispin)
                gj(k,l,ispin) = 0.0
             ELSE
                IF (l_flag) THEN
                   DO jspin = 1, jspins
                      fj(k,l,jspin) = ws(jspin) * ( uds(l,jspin)*gg - duds(l,jspin)*ff )
                      gj(k,l,jspin) = ws(jspin) * ( dus(l,jspin)*ff - us(l,jspin)*gg )
                   END DO
                ELSE
                   fj(k,l,ispin) = ws(ispin) * ( uds(l,ispin)*gg - duds(l,ispin)*ff )
                   gj(k,l,ispin) = ws(ispin) * ( dus(l,ispin)*ff - us(l,ispin)*gg )
                ENDIF
             ENDIF
          ENDDO
  ENDDO ! k = 1, lapw%nv

  END SUBROUTINE synth_fjgj


  SUBROUTINE hsmt_fjgj_gpu(input,atoms,cell,lapw,noco,usdus,n,ispin,fj,gj)
    !Calculate the fj&gj array which contain the part of the A,B matching coeff. depending on the
    !radial functions at the MT boundary as contained in usdus
    USE m_constants, ONLY : fpi_const
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
  
    REAL,INTENT(OUT),MANAGED     :: fj(:,0:,:,:),gj(:,0:,:,:)
    !     ..
    !     .. Local Scalars ..
    REAL con1

    INTEGER l,lo,intspin
    LOGICAL l_socfirst
    !     .. Local Arrays ..
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
    DO intspin=1,MERGE(2,1,noco%l_noco)

       CALL synth_fjgj(lapw%nv(intspin),ispin,input%jspins,atoms%lmax(n),atoms%lmaxd,apw,noco%l_constr.or.l_socfirst,&
            lapw%rk(:,intspin),atoms%rmt(n),con1,usdus%uds(:,n,:),usdus%dus(:,n,:),usdus%us(:,n,:),usdus%duds(:,n,:),&
            fj(:,0:,:,intspin),gj(:,0:,:,intspin))

    ENDDO
    RETURN
  END SUBROUTINE hsmt_fjgj_gpu
#endif

  SUBROUTINE hsmt_fjgj_cpu(input,atoms,cell,lapw,noco,usdus,n,ispin,fj,gj)
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
  
    REAL,INTENT(OUT)     :: fj(:,0:,:,:),gj(:,0:,:,:)
    !     ..
    !     .. Local Scalars ..
    REAL con1,ff,gg,gs

    INTEGER k,l,lo,intspin,jspin
    LOGICAL l_socfirst
    !     .. Local Arrays ..
    REAL ws(input%jspins)
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
    DO intspin=1,MERGE(2,1,noco%l_noco)
       !$OMP PARALLEL DO DEFAULT(NONE) &
       !$OMP PRIVATE(l,gs,fb,gb,ws,ff,gg,jspin)&
       !$OMP SHARED(lapw,atoms,con1,usdus,l_socfirst,noco,input)&
       !$OMP SHARED(fj,gj,intspin,n,ispin,apw)
       DO k = 1,lapw%nv(intspin)
          gs = lapw%rk(k,intspin)*atoms%rmt(n)
          CALL sphbes(atoms%lmax(n),gs, fb)
          CALL dsphbs(atoms%lmax(n),gs,fb, gb)
!          !$OMP SIMD PRIVATE(ws,ff,gg)
          DO l = 0,atoms%lmax(n)
             !---> set up wronskians for the matching conditions for each ntype
             DO jspin = 1, input%jspins
                ws(jspin) = con1/(usdus%uds(l,n,jspin)*usdus%dus(l,n,jspin)&
                            - usdus%us(l,n,jspin)*usdus%duds(l,n,jspin))
             END DO
             ff = fb(l)
             gg = lapw%rk(k,intspin)*gb(l)
             IF ( apw(l) ) THEN
                fj(k,l,ispin,intspin) = 1.0*con1 * ff / usdus%us(l,n,ispin)
                gj(k,l,ispin,intspin) = 0.0
             ELSE
                IF (noco%l_constr.or.l_socfirst) THEN
                   DO jspin = 1, input%jspins
                      fj(k,l,jspin,intspin) = ws(jspin) * ( usdus%uds(l,n,jspin)*gg - usdus%duds(l,n,jspin)*ff )
                      gj(k,l,jspin,intspin) = ws(jspin) * ( usdus%dus(l,n,jspin)*ff - usdus%us(l,n,jspin)*gg )
                   END DO
                ELSE
                   fj(k,l,ispin,intspin) = ws(ispin) * ( usdus%uds(l,n,ispin)*gg - usdus%duds(l,n,ispin)*ff )
                   gj(k,l,ispin,intspin) = ws(ispin) * ( usdus%dus(l,n,ispin)*ff - usdus%us(l,n,ispin)*gg )
                ENDIF
             ENDIF
          ENDDO
!          !$OMP END SIMD
       ENDDO ! k = 1, lapw%nv
       !$OMP END PARALLEL DO
    ENDDO
    RETURN
  END SUBROUTINE hsmt_fjgj_cpu
END MODULE m_hsmt_fjgj
