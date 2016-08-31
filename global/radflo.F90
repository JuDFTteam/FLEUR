!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_radflo
  USE m_juDFT
CONTAINS
  SUBROUTINE radflo(atoms, ntyp,jsp,ello,vr, f,g,mpi, usdus,&
       uuilon,duilon,ulouilopn,flo,lout_all)
    !
    !***********************************************************************
    ! generates the scalar relativistic wavefunctions (flo) needed for the
    ! local orbitals at atom type n for angular momentum l.
    ! the values of the function and the radial derivative on the sphere
    ! (ulos,dulos) boundaries, the overlap with the other radial functions
    ! (uulon,dulon) and between the different local orbitals (uloulopn) are
    ! also calculated.
    !
    ! ulos    : the value of the radial function of a local orbital
    !           at the muffin tin radius
    ! dulos   : the value of the radial derivative of the radial function
    !           function of a local orbital at the muffin tin radius
    ! uulon   : overlap integral between the radial functions of a local
    !           obital and the flapw radial function with the same l
    ! dulon   : overlap integral between the radial functions of a local
    !           obital and the energy derivative of the flapw radial
    !           function with the same l
    ! uloulopn: overlap integral between the radial functions of two
    !           different local orbitals
    !           (only needed if they have the same l)
    ! l_dulo  : whether we use a dot u as local orbital (see setlomap.F)
    !
    ! p.kurz jul. 1996 gb 2001
    !
    ! ulo_der : specifies the order of the energy derivative to be used
    !           (0, 1, 2, ... for primitive, first, second derivatives etc.)
    ! uuilon  : overlap integral between the radial functions of the
    !           integral (multiplied by ulo_der) of a local orbital and the
    !           flapw radial function with the same l
    ! duilon  : overlap integral between the radial functions of the
    !           integral of a local orbital and the energy derivative of the
    !           flapw radial function with the same l
    ! ulouilopn: overlap integral between the radial functions of the
    !           integral of a local orbital and another local orbital with
    !           the same l.
    !
    ! C. Friedrich Feb. 2005
    !***********************************************************************
    !
    USE m_intgr, ONLY : intgr0
    USE m_constants, ONLY : c_light
    USE m_radsra
    USE m_radsrdn
    USE m_differ
#include "cpp_double.h"
    USE m_types
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(INOUT):: usdus !lo part is calculated here
    TYPE(t_mpi),INTENT(IN)     :: mpi
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp,ntyp 
    LOGICAL, INTENT (IN), OPTIONAL :: lout_all
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: ello(atoms%nlod,atoms%ntypd),vr(atoms%jmtd)
    REAL,    INTENT (IN) :: f(atoms%jmtd,2,0:atoms%lmaxd),g(atoms%jmtd,2,0:atoms%lmaxd)
    REAL,    INTENT (OUT):: uuilon(atoms%nlod,atoms%ntypd),duilon(atoms%nlod,atoms%ntypd)
    REAL,    INTENT (OUT):: ulouilopn(atoms%nlod,atoms%nlod,atoms%ntypd)
    REAL,    INTENT (OUT):: flo(atoms%jmtd,2,atoms%nlod)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,j,k,l,ilo,jlo,nodelo,noded ,ierr,msh
    REAL    rn,t1,t2,d ,rr,fn,fl,fj,e,uds,duds,ddn,c,rmt
    LOGICAL ofdiag, loutput
    !     ..
    !     .. Local Arrays ..
    REAL dulo(atoms%jmtd),ulo(atoms%jmtd),glo(atoms%jmtd,2),filo(atoms%jmtd,2,atoms%nlod)
    REAL, ALLOCATABLE :: f_rel(:,:),vrd(:)
    REAL help(atoms%nlod+2,atoms%nlod+2)
    !     ..
    c = c_light(1.0)
    !
    IF ( PRESENT(lout_all) ) THEN
       loutput = ( mpi%irank == 0 ) .OR. lout_all
    ELSE
       loutput = ( mpi%irank == 0 )
    END IF
    !$    loutput=.false.
    IF (loutput) THEN
       WRITE (6,FMT=8000)
    END IF
    ofdiag = .FALSE.
    !---> calculate the radial wavefunction with the appropriate
    !---> energy parameter (ello)
    DO ilo = 1,atoms%nlo(ntyp)
       CALL radsra(ello(ilo,ntyp),atoms%llo(ilo,ntyp),vr(:),atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c,&
            usdus%ulos(ilo,ntyp,jsp),usdus%dulos(ilo,ntyp,jsp),nodelo,flo(:,1,ilo), flo(:,2,ilo))
       !
       !+apw+lo
       IF (atoms%l_dulo(ilo,ntyp).OR.atoms%ulo_der(ilo,ntyp).GE.1) THEN
          IF (atoms%ulo_der(ilo,ntyp).LE.8) THEN
             !--->    calculate orthogonal energy derivative at e
             i = atoms%ulo_der(ilo,ntyp)
             IF(atoms%l_dulo(ilo,ntyp)) i=1
             CALL radsrdn(ello(ilo,ntyp),atoms%llo(ilo,ntyp),vr(:),atoms%rmsh(1,ntyp), atoms%dx(ntyp),&
                  atoms%jri(ntyp),c, uds,duds,ddn,noded,glo,filo(:,:,ilo), flo(:,:,ilo),usdus%dulos(ilo,ntyp,jsp),i)
             DO i=1,atoms%jri(ntyp)
                flo(i,1,ilo) = glo(i,1)
                flo(i,2,ilo) = glo(i,2)
             ENDDO
             nodelo = noded
             ddn    = SQRT(ddn)
             IF(atoms%l_dulo(ilo,ntyp)) ddn=1.0
             flo (:,:,ilo)   = flo (:,:,ilo)/ddn ! Normalize ulo (flo) if APW+lo is not used
             filo(:,:,ilo)   = filo(:,:,ilo)/ddn ! and scale its integral (filo) accordingly
             usdus%dulos(ilo,ntyp,jsp) = duds/ddn          !   (setabc1lo and slomat assume <flo|flo>=1)
             usdus%ulos (ilo,ntyp,jsp) =  uds/ddn          ! 
          ELSE
             !
             !          test:
             !
             ! set up core-mesh
             d = EXP(atoms%dx(ntyp))
             rmt = atoms%rmsh(1,ntyp)
             DO i = 1, atoms%jri(ntyp) - 1
                rmt = rmt * d
             ENDDO
             rn = rmt
             msh = atoms%jri(ntyp)
             DO WHILE (rn < rmt + 20.0)
                msh = msh + 1
                rn = rn*d
             ENDDO
             rn = atoms%rmsh(1,ntyp)*( d**(msh-1) )
             ALLOCATE ( f_rel(msh,2),vrd(msh) )

             ! extend core potential (linear with slope t1 / a.u.)

             DO j = 1, atoms%jri(ntyp)
                vrd(j) = vr(j)
             ENDDO
             t1=0.125
             t2 = vrd(atoms%jri(ntyp))/rmt - rmt*t1
             rr = rmt
             DO j = atoms%jri(ntyp) + 1, msh
                rr = d*rr
                vrd(j) = rr*( t2 + rr*t1 )
             ENDDO
             e = ello(ilo,ntyp)
             fn = 6.0 ; fl = 1.0 ; fj = 0.5
             CALL differ(fn,fl,fj,c,82.0,atoms%dx(ntyp),atoms%rmsh(1,ntyp),&
                  rn,d,msh,vrd, e, f_rel(1,1),f_rel(1,2),ierr)

             f_rel(:,1) = 2 * f_rel(:,1)
             f_rel(:,2) = 2 * f_rel(:,2)
             rn = atoms%rmsh(1,ntyp)
             DO i = 1, atoms%jri(ntyp) 
                rn = rn * d
                WRITE(123,'(5f20.15)') rn,f_rel(i,1),f_rel(i,2), f(i,1,1),f(i,2,1)
             ENDDO

             STOP
          ENDIF

       ENDIF
       !-apw+lo
       !
       !--->    calculate the overlap between these fcn. and the radial functions
       !--->    of the flapw basis with the same l
       DO i = 1,atoms%jri(ntyp)
          ulo(i) = f(i,1,atoms%llo(ilo,ntyp))*flo(i,1,ilo) + f(i,2,atoms%llo(ilo,ntyp))*flo(i,2,ilo)
          dulo(i) = g(i,1,atoms%llo(ilo,ntyp))*flo(i,1,ilo) + g(i,2,atoms%llo(ilo,ntyp))*flo(i,2,ilo)
       END DO
       CALL intgr0(ulo, atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),usdus%uulon(ilo,ntyp,jsp))
       CALL intgr0(dulo,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),usdus%dulon(ilo,ntyp,jsp))
       IF (atoms%l_dulo(ilo,ntyp)) usdus%dulon(ilo,ntyp,jsp) = 0.0
       IF (loutput) THEN
          WRITE (6,FMT=8010) ilo,atoms%llo(ilo,ntyp),ello(ilo,ntyp),&
               usdus%ulos(ilo,ntyp,jsp),usdus%dulos(ilo,ntyp,jsp),nodelo,usdus%uulon(ilo,ntyp,jsp),&
               usdus%dulon(ilo,ntyp,jsp)
       END IF
       !
       !--->   case LO = energy derivative (ulo_der>=1):
       !--->   calculate the overlap between the LO-integral (filo) and the radial functions
       IF(atoms%ulo_der(ilo,ntyp).GE.1) THEN
          DO i=1,atoms%jri(ntyp)
             ulo(i) = f(i,1,atoms%llo(ilo,ntyp))*filo(i,1,ilo) + f(i,2,atoms%llo(ilo,ntyp))*filo(i,2,ilo)
             dulo(i) = g(i,1,atoms%llo(ilo,ntyp))*filo(i,1,ilo) + g(i,2,atoms%llo(ilo,ntyp))*filo(i,2,ilo)
          ENDDO
          CALL intgr0(ulo, atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),uuilon(ilo,ntyp))
          CALL intgr0(dulo,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),duilon(ilo,ntyp))
       ELSE
          uuilon(ilo,ntyp)=0
          duilon(ilo,ntyp)=0
       ENDIF
       !--->   calculate overlap between radial fcn. of different local
       !--->   orbitals (only if both have the same l)
       !         uloulopn(ilo,ilo,ntyp) = 1.0
       !         DO jlo = 1, (ilo-1)
       DO jlo = 1, ilo
          IF (atoms%llo(ilo,ntyp).EQ.atoms%llo(jlo,ntyp)) THEN
             DO i = 1,atoms%jri(ntyp)
                ulo(i) = flo(i,1,ilo)*flo(i,1,jlo) + flo(i,2,ilo)*flo(i,2,jlo)
             END DO
             CALL intgr0(ulo,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),usdus%uloulopn(ilo,jlo,ntyp,jsp))
             usdus%uloulopn(jlo,ilo,ntyp,jsp) = usdus%uloulopn(ilo,jlo,ntyp,jsp)
             ofdiag = .TRUE.
          ELSE
             usdus%uloulopn(ilo,jlo,ntyp,jsp) = 0.0
             usdus%uloulopn(jlo,ilo,ntyp,jsp) = 0.0
          END IF
       END DO
    END DO
    !
    !---> case: one of LOs = energy derivative (ulo_der>=1):
    !---> calculate overlap between LOs and integrals of LOs
    DO ilo = 1,atoms%nlo(ntyp)
       DO jlo = 1,atoms%nlo(ntyp)
          IF(atoms%ulo_der(jlo,ntyp).GE.1.AND.atoms%llo(ilo,ntyp).EQ.atoms%llo(jlo,ntyp)) THEN
             DO i = 1,atoms%jri(ntyp)
                ulo(i) = flo(i,1,ilo)*filo(i,1,jlo) + flo(i,2,ilo)*filo(i,2,jlo)
             ENDDO
             CALL intgr0(ulo,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),ulouilopn(ilo,jlo,ntyp))
          ELSE
             ulouilopn(ilo,jlo,ntyp)=0.0
          ENDIF
       ENDDO
    ENDDO
    ! 
    IF ( (ofdiag).AND.(loutput) ) THEN
       WRITE (6,FMT=*) 'overlap matrix between different local orbitals'
       WRITE (6,FMT=8020) (i,i=1,atoms%nlo(ntyp))
       DO ilo = 1,atoms%nlo(ntyp)
          WRITE (6,FMT='(i3,40e13.6)') ilo, (usdus%uloulopn(ilo,jlo,ntyp,jsp),jlo=1,atoms%nlo(ntyp))

       END DO
    END IF
    !
    !
    !     Diagonalize overlap matrix of normalized MT functions
    !       to check linear dependencies
    !       help is overlap matrix of all MT functions (flapw & LO)
    IF (.FALSE.) THEN
       DO l=0,MAXVAL(atoms%llo(1:atoms%nlo(ntyp),ntyp))
          IF(ALL(atoms%llo(1:atoms%nlo(ntyp),ntyp).NE.l)) CYCLE
          help(1,1)=1.0
          help(1,2)=0.0
          help(2,1)=0.0
          help(2,2)=1.0
          DO i=1,atoms%jri(ntyp)
             ulo(i)=g(i,1,l)**2+g(i,2,l)**2
          ENDDO
          CALL intgr0(ulo,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),ddn)
          ddn=SQRT(ddn)
          j=2
          k=2
          DO ilo=1,atoms%nlo(ntyp)
             IF(atoms%llo(ilo,ntyp).EQ.l) THEN
                j=j+1
                help(1,j)=usdus%uulon(ilo,ntyp,jsp)
                help(2,j)=usdus%dulon(ilo,ntyp,jsp)/ddn
                k=2
                DO jlo=1,ilo
                   IF(atoms%llo(jlo,ntyp).EQ.l) THEN
                      k=k+1
                      help(k,j)=usdus%uloulopn(ilo,jlo,ntyp,jsp)
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
          IF ( loutput ) THEN
             WRITE(6,'(A,I2)')&
                  &      'Overlap matrix of normalized MT functions '//&
                  &      'and its eigenvalues for l =',l
             DO i=1,j
                WRITE(6,'(30F11.7)') (help(k,i),k=1,i)
             ENDDO
          END IF
          CALL CPP_LAPACK_ssyev('N','U',j,help,atoms%nlod+2,ulo,dulo, i)
          IF(i/=0)  CALL juDFT_error("ssyev failed.",calledby ="radflo")
          IF ( loutput ) THEN
             WRITE(6,'(30F11.7)') (ulo(i),i=1,j)
          END IF
       ENDDO
    ENDIF ! .false.
    !

8000 FORMAT (/,t20,'radial function for local orbitals',/,t2,'lo',t6,&
         &       'l',t11,'energy',t29,'value',t42,'derivative',t56,'nodes',&
         &       t63,'ovlp with u',t78,'ovlp with udot')
8010 FORMAT (i3,i3,f10.5,5x,1p,2e16.7,i5,2e16.7)
8020 FORMAT (' lo',i9,19i15)
    RETURN
  END SUBROUTINE radflo
END MODULE m_radflo
