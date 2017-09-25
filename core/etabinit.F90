!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_etabinit
  USE m_juDFT
  !     *******************************************************
  !     *****   set up etab via old core program          *****
  !     *******************************************************
  !     modified to run with core-levels as provided by setcor
  !     ntab & ltab transport this info to core.F        gb`02
  !------------------------------------------------------------
CONTAINS
  SUBROUTINE etabinit(atoms,DIMENSION,input, vr,&
       etab,ntab,ltab,nkmust)

    USE m_constants, ONLY : c_light
    USE m_setcor
    USE m_differ
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_input),INTENT(IN)       :: input
    !
    !     .. Scalar Arguments ..
    !     ..
    !     .. Array Arguments ..
    REAL   , INTENT (IN) :: vr(atoms%jmtd,atoms%ntype) 
    REAL   , INTENT (OUT):: etab(100,atoms%ntype)
    INTEGER, INTENT (OUT):: ntab(100,atoms%ntype),ltab(100,atoms%ntype)
    INTEGER, INTENT (OUT):: nkmust(atoms%ntype)
    !     ..
    !     .. Local Scalars ..
    REAL  c,d,dxx,e,fj,fl,fn,rn,rnot,t2 ,z,t1,rr,weight
    REAL  bmu
    INTEGER i,ic,iksh,ilshell,j,jatom,korb,l, nst,ncmsh ,nshell,ipos,ierr
    !     ..
    !     .. Local Arrays ..
    INTEGER kappa(DIMENSION%nstd),nprnc(DIMENSION%nstd)
    REAL eig(DIMENSION%nstd),occ(DIMENSION%nstd,1),vrd(DIMENSION%msh),a(DIMENSION%msh),b(DIMENSION%msh)
    !     ..
    !
    c = c_light(1.0)
    !
    WRITE (6,FMT=8020)
    !
    ncmsh = DIMENSION%msh
    !     ---> set up densities
    DO  jatom = 1,atoms%ntype
       z = atoms%zatom(jatom)
       rn = atoms%rmt(jatom)
       dxx = atoms%dx(jatom)
       bmu = 0.0
       CALL setcor(jatom,1,atoms,input,bmu,nst,kappa,nprnc,occ)
       rnot = atoms%rmsh(1,jatom)
       d = EXP(atoms%dx(jatom))
       rn = rnot* (d** (ncmsh-1))
       WRITE (6,FMT=8000) z,rnot,dxx,atoms%jri(jatom)
       WRITE (16,FMT=8000) z,rnot,dxx,atoms%jri(jatom)
       DO j = 1,atoms%jri(jatom)
          vrd(j) = vr(j,jatom)
       ENDDO
       IF (input%l_core_confpot) THEN
          !--->    linear extension of the potential with slope t1 / a.u.
          t1=0.125
          t2  = vrd(atoms%jri(jatom))/atoms%rmt(jatom)-atoms%rmt(jatom)*t1
          rr = atoms%rmt(jatom)
          d = EXP(atoms%dx(jatom))
       ELSE
          t2 = vrd(atoms%jri(jatom))/ (atoms%jri(jatom)-DIMENSION%msh)
       ENDIF
       IF (atoms%jri(jatom).LT.DIMENSION%msh) THEN
          DO i = atoms%jri(jatom) + 1,DIMENSION%msh
             if (input%l_core_confpot) THEN
                rr = d*rr
                vrd(i) = rr*( t2 + rr*t1 )
             ELSE
                
                vrd(i) = vrd(atoms%jri(jatom)) + t2* (i-atoms%jri(jatom))
             ENDIF
          ENDDO
       END IF

       nst = atoms%ncst(jatom)
       DO  korb = 1,nst
          fn = nprnc(korb)
          fj = iabs(kappa(korb)) - .5e0
          weight = 2*fj + 1.e0
          fl = fj + (.5e0)*isign(1,kappa(korb))
          e = -2* (z/ (fn+fl))**2
          CALL differ(fn,fl,fj,c,z,dxx,rnot,rn,d,DIMENSION%msh,vrd,&
               e, a,b,ierr)
          IF (ierr/=0)  CALL juDFT_error("error in core-levels",calledby="etabinit")
          WRITE (6,FMT=8010) fn,fl,fj,e,weight
          WRITE (16,FMT=8010) fn,fl,fj,e,weight
          eig(korb) = e
       ENDDO
       ic = 0
       DO korb = 1,nst
          fn = nprnc(korb)
          fj = iabs(kappa(korb)) - .5e0
          weight = 2*fj + 1.e0
          fl = fj + (.5e0)*isign(1,kappa(korb))
          DO i = 1, INT(weight)
             ic = ic + 1
             IF (kappa(korb).GT.0) THEN
                ipos = ic + 1 + i 
             ELSEIF (kappa(korb).LT.-1) THEN
                ipos = ic - 2*(iabs(kappa(korb))-1) + MAX(i-2,0)
             ELSE
                ipos = ic
             ENDIF
             etab(ipos,jatom) = eig(korb)
             ntab(ipos,jatom) = NINT(fn)
             ltab(ipos,jatom) = NINT(fl)
          ENDDO
       ENDDO
       nkmust(jatom) = ic

       DO i=1,nkmust(jatom)
          WRITE(6,'(f12.6,2i3)') etab(i,jatom),ntab(i,jatom), ltab(i,jatom)
       ENDDO

    ENDDO
8000 FORMAT (/,/,10x,'z=',f4.0,5x,'r(1)=',e14.6,5x,'dx=',f8.6,5x,&
                'm.t.index=',i4,/,15x,'n',4x,'l',5x,'j',4x,'energy',7x, 'weight')
8010 FORMAT (12x,2f5.0,f6.1,f10.4,f12.4)
8020 FORMAT (/,/,12x,'core e.v. initialization')

  END SUBROUTINE etabinit
END MODULE m_etabinit
