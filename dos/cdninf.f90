!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdninf
CONTAINS
  SUBROUTINE cdninf(&
       &                  input,sym,noco,jspin,atoms,vacuum,&
       &                  sliceplot,banddos,ikpt,bkpt,wk,&
       &                  cell,kpts,&
       &                  nbands,eig,qal,qis,qvac,qvlay,&
       &                  qstars,jsym,ksym)
    !***********************************************************************
    !     this subroutine calculates the charge distribution of each state
    !     and writes this information to the out file. If dos or vacdos
    !     are .true. it also write the necessary information for dos or
    !     bandstructure plots to the file dosinp and vacdos respectivly
    !***********************************************************************
    !       changed this subroutine slightly for parallisation of dosinp&
    !       vacdos output (argument z replaced by ksym,jsym, removed sympsi
    !       call)                                        d.wortmann 5.99
    !
    !******** ABBREVIATIONS ************************************************
    !     qal      : l-like charge of each state
    !     qvac     : vacuum charge of each state
    !     qvlay    : charge in layers (z-ranges) in the vacuum of each state
    !     starcoeff: T if star coefficients have been calculated
    !     qstars   : star coefficients for layers (z-ranges) in vacuum
    !
    !***********************************************************************
    USE m_types
    IMPLICIT NONE
    TYPE(t_banddos),INTENT(IN)     :: banddos
    TYPE(t_sliceplot),INTENT(IN)   :: sliceplot
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms

    !     .. Scalar Arguments ..
    REAL,INTENT(IN):: wk
    INTEGER,INTENT(IN):: ikpt,jspin ,nbands 
    !
    !     STM Arguments
    COMPLEX, INTENT (IN) ::qstars(:,:,:,:) !(vacuum%nstars,DIMENSION%neigd,vacuum%layerd,2)
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: qvlay(:,:,:) !DIMENSION%neigd,vacuum%layerd,2)
    REAL,    INTENT (IN) :: qis(:,:,:)!(DIMENSION%neigd,kpts%nkptd,DIMENSION%jspd) 
    REAL,    INTENT (IN) :: qvac(:,:,:,:) !(DIMENSION%neigd,2,kpts%nkptd,DIMENSION%jspd)
    REAL,    INTENT (IN) :: bkpt(3)
    REAL,    INTENT (IN) :: eig(:)!(DIMENSION%neigd)
    REAL,    INTENT (IN) :: qal(0:,:,:)!(0:3,atoms%ntypd,neigd)
    INTEGER, INTENT (IN) :: jsym(:)!(DIMENSION%neigd)
    INTEGER, INTENT (IN) :: ksym(:)!(neigd)
    !     ..
    !     .. Local Scalars ..
    REAL qalmax,qishlp,qvacmt,qvact
    INTEGER i,iband,ilay,iqispc,iqvacpc,ityp,itypqmax,ivac,l,lqmax
    INTEGER istar
    !     ..
    !     .. Local Arrays ..
    REAL cartk(3)
    INTEGER iqalpc(0:3,atoms%ntypd)
    CHARACTER chstat(0:3)
    !     ..
    !     .. Data statements ..
    DATA chstat/'s','p','d','f'/
    !     ..


    IF (input%film) THEN
       WRITE (6,FMT=8000) (bkpt(i),i=1,3)
       WRITE (16,FMT=8000) (bkpt(i),i=1,3)
8000   FORMAT (/,3x,'q(atom,l): k=',3f10.5,/,/,t8,'e',t13,'max',t18,&
            &          'int',t22,'vac',t28,'spheres(s,p,d,f)')
       IF (banddos%dos) THEN
          cartk=MATMUL(bkpt,cell%bmat)
          WRITE (85,FMT=8020) cartk(1),cartk(2),cartk(3),nbands,wk
          !     *************** for vacdos shz Jan.96
          IF (banddos%vacdos) THEN
             WRITE (86,FMT=8020) cartk(1),cartk(2),cartk(3),nbands,wk
          END IF
       END IF
    ELSE
       WRITE (6,FMT=8010) (bkpt(i),i=1,3)
       WRITE (16,FMT=8010) (bkpt(i),i=1,3)
8010   FORMAT (/,3x,'q(atom,l): k=',3f10.5,/,/,t8,'e',t13,'max',t18,&
            &          'int',t24,'spheres(s,p,d,f)')
       IF (banddos%dos) THEN
          cartk=MATMUL(bkpt,cell%bmat)
          WRITE (85,FMT=8020) cartk(1),cartk(2),cartk(3),nbands,wk
       END IF
    END IF
8020 FORMAT (1x,3e20.12,i6,e20.12)

    DO iband = 1,nbands
       IF (sliceplot%slice) THEN
          WRITE (6,FMT=8030) iband,eig(iband)
          WRITE (16,FMT=8030) iband,eig(iband)
8030      FORMAT (' cdnval: slice for i=',i4,'  and energy',1e12.4)
       END IF

       qvacmt = 0.0
       qvact = 0.0
       IF (input%film) THEN
          DO ivac = 1,vacuum%nvac
             qvact = qvact + qvac(iband,ivac,ikpt,jspin)
          END DO
          IF (sym%invs .OR. sym%zrfs) qvact = 2.0*qvact
          iqvacpc = NINT(qvact*100.0)
          qvacmt = qvact
       END IF
       qalmax = 0.0
       lqmax = 0
       itypqmax = 0
       DO ityp = 1,atoms%ntype
          DO l = 0,3
             iqalpc(l,ityp) = NINT(qal(l,ityp,iband)*100.0)
             qvacmt = qvacmt + qal(l,ityp,iband)*atoms%neq(ityp)
             IF (qalmax.LT.qal(l,ityp,iband)) THEN
                qalmax = qal(l,ityp,iband)
                lqmax = l
                itypqmax = ityp
             END IF
          END DO
       END DO
       qishlp = 1.0 - qvacmt
       IF (noco%l_noco) qishlp = qis(iband,ikpt,jspin)
       iqispc = NINT(qishlp*100.0)
       IF (input%film) THEN
          WRITE (6,FMT=8040) eig(iband),chstat(lqmax),itypqmax,&
               &        iqispc,iqvacpc, ((iqalpc(l,ityp),l=0,3),ityp=1,atoms%ntype)
          WRITE (16,FMT=8040) eig(iband),chstat(lqmax),itypqmax,&
               &        iqispc,iqvacpc, ((iqalpc(l,ityp),l=0,3),ityp=1,atoms%ntype)
8040      FORMAT (f10.4,2x,a1,i2,2x,2i3, (t26,6 (4i3,1x)))
          IF (banddos%dos) THEN
             IF (banddos%ndir.NE.0) THEN
                WRITE (85,FMT=8050) eig(iband),ksym(iband),&
                     &              jsym(iband),qvact, ((qal(l,ityp,iband),l=0,3),&
                     &              ityp=1,atoms%ntype)
8050            FORMAT (f12.5,2i2,f12.5,/, (4f12.5))
             ELSE
                WRITE (85,FMT=8060) eig(iband),&
                     &              ((qal(l,ityp,iband),l=0,3),ityp=1,atoms%ntype),qvact
8060            FORMAT (10f12.7)
             END IF
          END IF
          !     ***************** for vacdos shz Jan.96
          IF (banddos%vacdos) THEN
             IF (.NOT.vacuum%starcoeff) THEN
                WRITE (86,FMT=8070) eig(iband),&
                     &               ((qvlay(iband,ilay,ivac),ilay=1,vacuum%layers),&
                     &                               ivac=1,vacuum%nvac)
             ELSE
                WRITE (86,FMT=8070) eig(iband),&
                     &                 ((qvlay(iband,ilay,ivac),&
                     &                 (REAL(qstars(istar,iband,ilay,ivac)),&
                     &                 istar=1,vacuum%nstars-1),ilay=1,vacuum%layers),ivac=1,vacuum%nvac)
             END IF
8070         FORMAT (f10.4,2x,20(e16.8,1x))

          END IF
          !     **************************************
       ELSE
          WRITE (6,FMT=8080) eig(iband),chstat(lqmax),itypqmax,&
               &        iqispc, ((iqalpc(l,ityp),l=0,3),ityp=1,atoms%ntype)
          WRITE (16,FMT=8080) eig(iband),chstat(lqmax),itypqmax,&
               &        iqispc, ((iqalpc(l,ityp),l=0,3),ityp=1,atoms%ntype)
8080      FORMAT (f10.4,2x,a1,i2,2x,i3, (t26,6 (4i3,1x)))
          IF (banddos%dos) THEN
             IF (banddos%ndir.NE.0) THEN
                WRITE (85,FMT=8050) eig(iband),ksym(iband),&
                     &              jsym(iband),0.0, ((qal(l,ityp,iband),l=0,3),ityp=1,&
                     &              atoms%ntype)
             ELSE
                WRITE (85,FMT=8060) eig(iband),&
                     &              ((qal(l,ityp,iband),l=0,3),ityp=1,atoms%ntype),0.0
             END IF
          END IF
       END IF
    END DO
  END SUBROUTINE cdninf
END MODULE m_cdninf
