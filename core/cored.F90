MODULE m_cored
CONTAINS
  SUBROUTINE cored(&
       &                 input,jspin,atoms,&
       &                 rho,DIMENSION,&
       &                 sphhar,&
       &                 vr,&
       &                 qint,rhc,seig)

    !     *******************************************************
    !     *****   set up the core densities for compounds.  *****
    !     *****                      d.d.koelling           *****
    !     *******************************************************
    USE m_juDFT
    USE m_intgr, ONLY : intgr3,intgr0,intgr1
    USE m_constants, ONLY : c_light,sfp_const
    USE m_setcor
    USE m_differ
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
    !
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jspin   
    REAL,    INTENT (OUT) :: seig
    !     ..
    !     .. Array Arguments ..
    REAL   , INTENT (IN) :: vr(atoms%jmtd,atoms%ntypd)
    REAL,    INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,DIMENSION%jspd)
    REAL,    INTENT (OUT) :: rhc(DIMENSION%msh,atoms%ntypd),qint(atoms%ntypd,DIMENSION%jspd)
    !     ..
    !     .. Local Scalars ..
    REAL e,fj,fl,fn,q,rad,rhos,rhs,sea,sume,t2,tec 
    REAL d,dxx,rn,rnot,z,t1,rr,r,lambd,c,bmu,weight
    INTEGER i,j,jatom,jm,korb,n,ncmsh,nm,nm1,nst ,l,ierr
    !     ..
    !     .. Local Arrays ..
    REAL rhcs(DIMENSION%msh),rhoc(DIMENSION%msh),rhoss(DIMENSION%msh),vrd(DIMENSION%msh),f(0:3)
    REAL occ(DIMENSION%nstd),a(DIMENSION%msh),b(DIMENSION%msh),ain(DIMENSION%msh),ahelp(DIMENSION%msh)
    REAL occ_h(DIMENSION%nstd,2)
    INTEGER kappa(DIMENSION%nstd),nprnc(DIMENSION%nstd)
    !     ..
    c = c_light(1.0)
    seig = 0.
    IF (jspin.EQ.1) THEN
       OPEN (17,file='cdnc',form='unformatted',status='unknown')
    ENDIF
    !
    IF (input%frcor) THEN
       IF (jspin.EQ.1) REWIND 17
       DO  n = 1,atoms%ntype
          jm = atoms%jri(n)
          rnot = atoms%rmsh(1,n) ; dxx = atoms%dx(n)
          ncmsh = NINT( LOG( (atoms%rmt(n)+10.0)/rnot ) / dxx + 1 )
          ncmsh = MIN( ncmsh, DIMENSION%msh )
          !     --->    read in core density
          READ (17) (rhc(i,n),i=1,ncmsh)
          !     --->    update spherical charge density
          DO  i = 1,atoms%jri(n)
             rhoc(i) = rhc(i,n)
             rho(i,0,n,jspin) = rho(i,0,n,jspin) + rhoc(i)/sfp_const
          ENDDO
          !     --->    read in kinetic enrgy of the core
          READ (17) tec
          !     ---> for total energy calculations, determine the sum of the
          !     ---> eigenvalues by requiring that the core kinetic energy
          !     ---> remains constant.
          DO  i = 1,atoms%jri(n)
             rhoc(i) = rhoc(i)*vr(i,n)/atoms%rmsh(i,n)
          ENDDO
          nm = atoms%jri(n)
          CALL intgr3(rhoc,atoms%rmsh(1,n),atoms%dx(n),nm,rhos)
          sea = tec + rhos
          WRITE (16,FMT=8030) n,jspin,tec,sea
          WRITE (6,FMT=8030) n,jspin,tec,sea
          seig = seig + atoms%neq(n)*sea
       ENDDO
       !     --->    read in qint
       READ (17) (qint(n,jspin),n=1,atoms%ntype)
       RETURN
    END IF

    !#ifdef CPP_CORE
    !      IF (jspin.EQ.1) THEN
    !        OPEN (45,file='slaterf',form='formatted',status='unknown')
    !      ENDIF
    !#endif
    !     ---> set up densities
    DO  jatom = 1,atoms%ntype
       sume = 0.
       z = atoms%zatom(jatom)
       !         rn = rmt(jatom)
       dxx = atoms%dx(jatom)
       bmu = 0.0
       CALL setcor(jatom,DIMENSION%jspd,atoms,bmu,nst,kappa,nprnc,occ_h)
       IF ((bmu > 99.)) THEN
          occ(1:nst) = input%jspins *  occ_h(1:nst,jspin)
       ELSE
          occ(1:nst) = occ_h(1:nst,1)
       ENDIF
       rnot = atoms%rmsh(1,jatom)
       d = EXP(atoms%dx(jatom))
       ncmsh = NINT( LOG( (atoms%rmt(jatom)+10.0)/rnot ) / dxx + 1 )
       ncmsh = MIN( ncmsh, DIMENSION%msh )
       rn = rnot* (d** (ncmsh-1))
       WRITE (6,FMT=8000) z,rnot,dxx,atoms%jri(jatom)
       WRITE (16,FMT=8000) z,rnot,dxx,atoms%jri(jatom)
       DO  j = 1,atoms%jri(jatom)
          rhoss(j) = 0.
          vrd(j) = vr(j,jatom)
       ENDDO
       !
#ifdef CPP_CORE
       !--->    linear extension of the potential with slope t1 / a.u.
       t1=0.125
       t1 = MAX( (vrd(atoms%jri(jatom)) - vrd(atoms%jri(jatom)-1)*d)*&
            d / (atoms%rmt(jatom)**2 * (d-1) ) , t1)
       t2=vrd(atoms%jri(jatom))/atoms%rmt(jatom)-atoms%rmt(jatom)*t1
       rr = atoms%rmt(jatom)
#else
       t2 = vrd(atoms%jri(jatom)) / ( atoms%jri(jatom) - ncmsh )
#endif
       IF ( atoms%jri(jatom) .LT. ncmsh) THEN
          DO  i = atoms%jri(jatom) + 1,ncmsh
             rhoss(i) = 0.
#ifdef CPP_CORE
             rr = d*rr
             vrd(i) = rr*( t2 + rr*t1 )
             !               vrd(i) = 2*vrd(jri(jatom)) - rr*( t2 + rr*t1 )
#else
             vrd(i) = vrd(atoms%jri(jatom)) + t2* (i-atoms%jri(jatom))
#endif
             !
          ENDDO
       END IF

       !#ifndef CPP_CORE
       nst = atoms%ncst(jatom)        ! for lda+U
       !#endif
       IF (input%gw.EQ.1 .OR. input%gw.EQ.3)&
            &                      WRITE(15) nst,atoms%rmsh(1:atoms%jri(jatom),jatom)

       DO  korb = 1,nst
          !#ifndef CPP_CORE
          IF (occ(korb).EQ.0) CYCLE
          !#endif
          fn = nprnc(korb)
          fj = iabs(kappa(korb)) - .5e0
          weight = 2*fj + 1.e0
          IF (bmu > 99.) weight = occ(korb)
          fl = fj + (.5e0)*isign(1,kappa(korb))
          e = -2* (z/ (fn+fl))**2
          CALL differ(fn,fl,fj,c,z,dxx,rnot,rn,d,ncmsh,vrd, e, a,b,ierr)
          WRITE (6,FMT=8010) fn,fl,fj,e,weight
          WRITE (16,FMT=8010) fn,fl,fj,e,weight
          IF (ierr/=0)  CALL juDFT_error("error in core-level routine" ,calledby ="cored")
          IF (input%gw.EQ.1 .OR. input%gw.EQ.3) WRITE (15) NINT(fl),weight,e,&
               a(1:atoms%jri(jatom)),b(1:atoms%jri(jatom))

          sume = sume + weight*e/input%jspins
          DO j = 1,ncmsh
             rhcs(j) = weight* (a(j)**2+b(j)**2)
             rhoss(j) = rhoss(j) + rhcs(j)
          ENDDO
          !#ifdef CPP_CORE
          !            ENDIF
          !#endif
       ENDDO

       !     ---->update spherical charge density rho with the core density.
       !     ---->for spin-polarized (jspins=2), take only half the density
       nm = atoms%jri(jatom)
       DO  j = 1,nm
          rhoc(j) = rhoss(j)/input%jspins
          rho(j,0,jatom,jspin) = rho(j,0,jatom,jspin) + rhoc(j)/sfp_const
       ENDDO

       rhc(1:ncmsh,jatom)   = rhoss(1:ncmsh) / input%jspins
       rhc(ncmsh+1:DIMENSION%msh,jatom) = 0.0

       seig = seig + atoms%neq(jatom)*sume
       !         WRITE (17) (rhoc(i),i=1,nm)
       WRITE (17) (rhc(i,jatom),i=1,ncmsh)
       DO  i = 1,nm
          rhoc(i) = rhoc(i)*vr(i,jatom)/atoms%rmsh(i,jatom)
       ENDDO
       CALL intgr3(rhoc,atoms%rmsh(1,jatom),atoms%dx(jatom),nm,rhs)
       tec = sume - rhs
       WRITE (6,FMT=8030) jatom,jspin,tec,sume
       WRITE (16,FMT=8030) jatom,jspin,tec,sume
       WRITE (17) tec

       !     ---> simpson integration
       rad = atoms%rmt(jatom)
       q = rad*rhoss(nm)/2.
       DO  nm1 = nm + 1,ncmsh - 1,2
          rad = d*rad
          q = q + 2*rad*rhoss(nm1)
          rad = d*rad
          q = q + rad*rhoss(nm1+1)
       ENDDO
       q = 2*q*dxx/3
       !+sb
       WRITE (6,FMT=8020) q/input%jspins
       WRITE (16,FMT=8020) q/input%jspins
       !-sb
       qint(jatom,jspin) = q*atoms%neq(jatom)

    ENDDO

#ifdef CPP_CORE
    IF (jspin.EQ.input%jspins) THEN
       CLOSE (45)
    ENDIF
#endif

    !      qint=0.
    WRITE (17) (qint(n,jspin),n=1,atoms%ntype)
    !
    IF (jspin.EQ.input%jspins) CLOSE (17)
    RETURN

8000 FORMAT (/,/,10x,'z=',f4.0,5x,'r(1)=',e14.6,5x,'dx=',f8.6,5x,&
         &       'm.t.index=',i4,/,15x,'n',4x,'l',5x,'j',4x,'energy',7x,&
         &       'weight')
8010 FORMAT (12x,2f5.0,f6.1,f10.4,f10.0)
8020 FORMAT (f20.8,'  electrons lost from core.')
8030 FORMAT (10x,'atom type',i3,'  (spin',i2,') ',/,10x,&
         &       'kinetic energy=',e20.12,5x,'sum of the eigenvalues=',&
         &       e20.12)
  END subroutine cored
END MODULE m_cored

