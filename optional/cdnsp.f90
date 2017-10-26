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
        SUBROUTINE cdnsp(&
             &                 atoms,input,vacuum,sphhar,&
             &                 stars,sym,oneD,cell,DIMENSION)

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
          TYPE(t_input),INTENT(INOUT)  :: input
          TYPE(t_sym),INTENT(IN)       :: sym
          TYPE(t_oneD),INTENT(IN)      :: oneD
          TYPE(t_cell),INTENT(IN)      :: cell
          TYPE(t_dimension),INTENT(IN) :: DIMENSION

          ! local type instances
          TYPE(t_potden)               :: den

          !     .. Local Scalars ..
          REAL dummy,p,pp,qtot1,qtot2,spmtot,qval,sfp,fermiEnergyTemp
          INTEGER i,ivac,j,k,lh,n,na,jsp_new
          INTEGER ios 
          LOGICAL n_exist,l_qfix
          !     ..
          !     .. Local Arrays ..
          REAL rhoc(atoms%jmtd,atoms%ntype,dimension%jspd)
          REAL tec(atoms%ntype,dimension%jspd),qintc(atoms%ntype,dimension%jspd)
          CHARACTER(len=140), ALLOCATABLE :: clines(:)
          CHARACTER(len=140)              :: lineread
          !      ..
          sfp = 2 * SQRT(pi_const)
          !sphhar%nlhd = MAXVAL(sphhar%nlh(:))

          IF (input%jspins/=2)  CALL juDFT_error&
               &     ("cdnsp: set jspins = 2 and remove fl7para!",calledby&
               &     ="cdnsp")

          CALL den%init(stars,atoms,sphhar,vacuum,oneD,DIMENSION%jspd,.FALSE.,POTDEN_TYPE_DEN)
          ALLOCATE (den%cdom(1),den%cdomvz(1,1),den%cdomvxy(1,1,1))
          ALLOCATE (den%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_u),input%jspins))
          den%mmpMat = CMPLX(0.0,0.0)

          input%jspins=1
          CALL readCoreDensity(input,atoms,dimension,rhoc,tec,qintc)
          CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,&
                           CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den%iter,den%mt,den%pw,den%vacz,&
                           den%vacxy,den%cdom,den%cdomvz,den%cdomvxy)
          input%jspins=2

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
             !         WRITE (16,FMT='(8f10.4)') (den%mt(i,0,n,1),i=1,16)
             CALL intgr3(den%mt(1,0,n,1),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),qval)
             p = (atoms%bmu(n)+sfp*qval)/ (2.*sfp*qval)
             pp = 1. - p
             DO j = 1,atoms%jri(n)
                den%mt(j,0,n,jsp_new) = pp*den%mt(j,0,n,1) + rhoc(j,n,1)/ (2.*sfp)
                den%mt(j,0,n,1)       =  p*den%mt(j,0,n,1) + rhoc(j,n,1)/ (2.*sfp)
             ENDDO
             DO lh = 1,sphhar%nlh(atoms%ntypsy(na))
                DO j = 1,atoms%jri(n)
                   den%mt(j,lh,n,jsp_new) = pp*den%mt(j,lh,n,1)
                   den%mt(j,lh,n,1)       =  p*den%mt(j,lh,n,1)
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
          !     ----> write the spin-polarized density
          CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,&
                            CDN_INPUT_DEN_const,0,-1.0,0.0,.FALSE.,den%iter,den%mt,den%pw,den%vacz,den%vacxy,&
                            den%cdom,den%cdomvz,den%cdomvxy)
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
             WRITE (6,'('' moment in sphere '',2x,'':'',f8.4)') spmtot
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
          !
          ! for lda+U: flip n-matrix
          !
          IF (atoms%n_u.GT.0) THEN
             INQUIRE (file='n_mmp_mat',exist=n_exist)
             IF (n_exist) THEN
                OPEN (69,file='n_mmp_mat',status='old',form='formatted')
                REWIND 69

                n=0
                DO
                   READ (69,'(a)',iostat=ios) lineread
                   IF (ios.NE.0) EXIT
                   n = n+1
                ENDDO
                ALLOCATE (clines(n))
                REWIND 69
                DO i=1,n
                   WRITE (69,'(a)') TRIM(clines(i))
                ENDDO
                DO i=1,n
                   WRITE (69,'(a)') TRIM(clines(i))
                ENDDO
                DEALLOCATE (clines)

                CLOSE(69)
             ENDIF
          ENDIF


        END SUBROUTINE cdnsp
      END MODULE m_cdnsp
