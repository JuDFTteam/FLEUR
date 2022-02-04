!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_cdnsp
      USE m_juDFT
!     *******************************************************
!     sets up the starting density for the spin-polarized
!     calculation from a paramagnetic density
!     changed to suit both ferromagnetic and antiferro-
!     magnetic case. changes only in mt-part - r.pentcheva Jan'96
!     *******************************************************
      CONTAINS
        SUBROUTINE cdnsp(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell)

          USE m_intgr, ONLY : intgr3
          USE m_constants
          USE m_cdn_io
          USE m_types
          IMPLICIT NONE
          !     ..
          TYPE(t_stars),INTENT(IN)     :: stars
          TYPE(t_vacuum),INTENT(IN)    :: vacuum
          TYPE(t_atoms),INTENT(IN)     :: atoms
          TYPE(t_sphhar),INTENT(IN)    :: sphhar
          TYPE(t_input),INTENT(IN)     :: input
          TYPE(t_sym),INTENT(IN)       :: sym
          TYPE(t_noco),INTENT(IN)      :: noco
          TYPE(t_oneD),INTENT(IN)      :: oneD
          TYPE(t_cell),INTENT(IN)      :: cell


          ! local type instances
          TYPE(t_potden)               :: den
          TYPE(t_input)                ::input_jsp
          !     .. Local Scalars ..
          REAL dummy,pp,qtot1,qtot2,spmtot,qval,sfp,fermiEnergyTemp,tempDistance
          INTEGER i,ivac,j,k,lh,n,na,jsp_new,i_u
          INTEGER ios, archiveType
          LOGICAL n_exist,l_qfix
          !     ..
          !     .. Local Arrays ..
          REAL p(atoms%ntype)
          REAL rhoc(atoms%jmtd,atoms%ntype,input%jspins)
          REAL tec(atoms%ntype,input%jspins),qintc(atoms%ntype,input%jspins)
          CHARACTER(len=140), ALLOCATABLE :: clines(:)
          CHARACTER(len=140)              :: lineread
          !      ..
          sfp = 2 * SQRT(pi_const)
          !sphhar%nlhd = MAXVAL(sphhar%nlh(:))

          IF (input%jspins/=2) CALL juDFT_error("cdnsp: set jspins = 2!", calledby ="cdnsp")

          CALL den%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
          input_jsp=input
          input_jsp%jspins=1
          CALL readCoreDensity(input_jsp,atoms,rhoc,tec,qintc)

          CALL readDensity(stars,noco,vacuum,atoms,cell,sphhar,input_jsp,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,&
                           CDN_INPUT_DEN_const,0,fermiEnergyTemp,tempDistance,l_qfix,den)

          qval = 0.
          na = 1
          !
          !     ---> set jspins=2
          jsp_new = 2
          !
          DO n = 1,atoms%ntype
             DO j = 1,atoms%jri(n)
                den%mt(j,0,n,1) = den%mt(j,0,n,1) - rhoc(j,n,1)/sfp
             ENDDO
             CALL intgr3(den%mt(1,0,n,1),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),qval)
             p(n) = (atoms%bmu(n)+sfp*qval)/ (2.*sfp*qval)
             pp = 1.0 - p(n)
             DO j = 1,atoms%jri(n)
                den%mt(j,0,n,jsp_new) = pp*den%mt(j,0,n,1) + rhoc(j,n,1)/ (2.*sfp)
                den%mt(j,0,n,1)       =  p(n)*den%mt(j,0,n,1) + rhoc(j,n,1)/ (2.*sfp)
             ENDDO
             DO lh = 1,sphhar%nlh(sym%ntypsy(na))
                DO j = 1,atoms%jri(n)
                   den%mt(j,lh,n,jsp_new) = pp*den%mt(j,lh,n,1)
                   den%mt(j,lh,n,1)       =  p(n)*den%mt(j,lh,n,1)
                ENDDO
             ENDDO
             na = na + atoms%neq(n)
          ENDDO
          DO k = 1,stars%ng3
             den%pw(k,jsp_new) = 0.5 * den%pw(k,1)
             den%pw(k,1)       = den%pw(k,jsp_new)
          ENDDO
          IF (input%film) THEN
             DO ivac = 1,vacuum%nvac
                DO j = 1, vacuum%nmz
                   den%vacz(j,ivac,jsp_new) = 0.5 * den%vacz(j,ivac,1)
                   den%vacz(j,ivac,1)       = den%vacz(j,ivac,jsp_new)
                ENDDO
                DO k = 2, stars%ng2
                   DO j = 1,vacuum%nmzxy
                      den%vacxy(j,k-1,ivac,jsp_new) = 0.5 * den%vacxy(j,k-1,ivac,1)
                      den%vacxy(j,k-1,ivac,1)       = den%vacxy(j,k-1,ivac,jsp_new)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          ! LDA + U
          IF (atoms%n_u.GT.0) THEN
             DO i_u = 1, atoms%n_u
                n = atoms%lda_u(i_u)%atomType
                pp = 1.0 - p(n)
                den%mmpMat(:,:,i_u,jsp_new) = pp * den%mmpMat(:,:,i_u,1)
                den%mmpMat(:,:,i_u,1) = p(n) * den%mmpMat(:,:,i_u,1)
             END DO
          END IF

          rhoc(:,:,1) = 0.5 * rhoc(:,:,1)
          rhoc(:,:,jsp_new) = rhoc(:,:,1)
          tec(:,1) = 0.5 * tec(:,1)
          tec(:,jsp_new) = tec(:,1)
          qintc(:,1) = 0.5 * qintc(:,1)
          qintc(:,jsp_new) = 0.5 * qintc(:,1)

          CALL writeCoreDensity(input,atoms,rhoc,tec,qintc)

          !     ----> write the spin-polarized density
          CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,&
                            CDN_INPUT_DEN_const,0,-1.0,0.0,-1.0,-1.0,.FALSE.,den)
          !
          !     -----> This part is only used for testing th e magnetic moment in
          !     ----->   each sphere
          !
          DO n = 1,atoms%ntype
             qtot1=0.00
             qtot2=0.00
             CALL intgr3(den%mt(1,0,n,1),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),qtot1)
             CALL intgr3(den%mt(1,0,n,jsp_new),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),qtot2)
             spmtot=sfp*(qtot1-qtot2)
             WRITE (oUnit,'('' moment in sphere '',2x,'':'',f8.4)') spmtot
          ENDDO

          !--->   read enpara and then double it
          INQUIRE(file='enpara',exist=n_exist)
          IF (n_exist) THEN
             OPEN(40,file ='enpara',status='old',form='formatted')
             REWIND 40
             n = 0
             DO
                READ (40,'(a)',iostat = ios) lineread
                IF (ios/=0) EXIT
                n          = n+1
             ENDDO

             ALLOCATE (clines(n))

             REWIND 40
             DO i = 1,n
                READ (40,'(a)') clines(i)
             ENDDO

             REWIND 40
             DO i = 1,n
                WRITE (40,'(a)') TRIM(clines(i))
             ENDDO
             DO i = 1,n
                WRITE (40,'(a)') TRIM(clines(i))
             ENDDO

             DEALLOCATE (clines)
             CLOSE(40)
          ENDIF
!          !
!          ! for lda+U: flip n-matrix
!          !
!          IF (atoms%n_u.GT.0) THEN
!             INQUIRE (file='n_mmp_mat',exist=n_exist)
!             IF (n_exist) THEN
!                OPEN (69,file='n_mmp_mat',status='old',form='formatted')
!                REWIND 69
!
!                n=0
!                DO
!                   READ (69,'(a)',iostat=ios) lineread
!                   IF (ios.NE.0) EXIT
!                   n = n+1
!                ENDDO
!                ALLOCATE (clines(n))
!                REWIND 69
!                DO i=1,n
!                   WRITE (69,'(a)') TRIM(clines(i))
!                ENDDO
!                DO i=1,n
!                   WRITE (69,'(a)') TRIM(clines(i))
!                ENDDO
!                DEALLOCATE (clines)
!
!                CLOSE(69)
!             ENDIF
!          ENDIF


        END SUBROUTINE cdnsp
      END MODULE m_cdnsp
