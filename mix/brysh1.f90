MODULE m_brysh1
  USE m_juDFT
  !******************************************************
  !      shifts the charge density of the interstitial, m.t
  !      and vacuum part in one single vector
  !      in the spin polarized case the arrays consist of 
  !      spin up and spin down densities
  !******************************************************
CONTAINS
  SUBROUTINE brysh1(input,stars,atoms,sphhar,noco,vacuum,sym,oneD,&
                    intfac,vacfac,den,nmap,nmaph,mapmt,mapvac,mapvac2,sout) 

    USE m_types
    IMPLICIT NONE
    TYPE(t_oneD),INTENT(IN)    :: oneD
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_potden),INTENT(IN)  :: den

    ! Scalar Arguments
    REAL,    INTENT (IN)  :: intfac,vacfac
    INTEGER, INTENT (OUT) :: mapmt,mapvac,mapvac2,nmap,nmaph

    ! Array Arguments
    REAL,    INTENT (OUT) :: sout(:)

    ! Local Scalars
    INTEGER i,iv,j,js,k,l,n,nall,na,nvaccoeff,nvaccoeff2,mapmtd

    !--->  put input into arrays sout 
    !      in the spin polarized case the arrays consist of 
    !      spin up and spin down densities

    j=0
    DO  js = 1,input%jspins
       DO i = 1,stars%ng3
          j = j + 1
          sout(j) = REAL(den%pw(i,js))
       END DO
       IF (.NOT.sym%invs) THEN
          DO i = 1,stars%ng3
             j = j + 1
             sout(j) = AIMAG(den%pw(i,js))
          END DO
       ENDIF
       mapmt=0
       na = 1
       DO n = 1,atoms%ntype
          DO l = 0,sphhar%nlh(atoms%ntypsy(na))
             DO i = 1,atoms%jri(n)
                mapmt = mapmt +1
                j = j + 1
                sout(j) = den%mt(i,l,n,js)
             END DO
          END DO
          na = na + atoms%neq(n)
       END DO
       mapvac=0
       IF (input%film) THEN
          DO iv = 1,vacuum%nvac
             DO k = 1,vacuum%nmz
                mapvac = mapvac + 1
                j = j + 1
                sout(j) = den%vacz(k,iv,js)
             END DO
             DO k = 1,oneD%odi%nq2-1
                DO i = 1,vacuum%nmzxy
                   mapvac = mapvac + 1
                   j = j + 1
                   sout(j) =  REAL(den%vacxy(i,k,iv,js))
                END DO
             END DO
             IF (.NOT.sym%invs2) THEN
                DO k = 1,oneD%odi%nq2-1
                   DO i = 1,vacuum%nmzxy
                      mapvac = mapvac + 1
                      j = j + 1
                      sout(j) =  AIMAG(den%vacxy(i,k,iv,js))
                   END DO
                END DO
             END IF
          END DO
       END IF
       IF (js .EQ. 1) nmaph = j
    ENDDO

    mapvac2=0
    IF (noco%l_noco) THEN
       !--->    off-diagonal part of the density matrix
       DO i = 1,stars%ng3
          j = j + 1
          sout(j) = REAL(den%cdom(i))
       END DO
       DO i = 1,stars%ng3
          j = j + 1
          sout(j) = AIMAG(den%cdom(i))
       END DO
       IF (input%film) THEN
          DO iv = 1,vacuum%nvac
             DO k = 1,vacuum%nmz
                mapvac2 = mapvac2 + 1
                j = j + 1
                sout(j) = REAL(den%cdomvz(k,iv))
             END DO
             DO k = 1,oneD%odi%nq2-1
                DO i = 1,vacuum%nmzxy
                   mapvac2 = mapvac2 + 1
                   j = j + 1
                   sout(j) =  REAL(den%vacxy(i,k,iv,3))
                END DO
             END DO
          END DO
          DO iv = 1,vacuum%nvac
             DO k = 1,vacuum%nmz
                mapvac2 = mapvac2 + 1
                j = j + 1
                sout(j) = AIMAG(den%cdomvz(k,iv))
             END DO
             DO k = 1,oneD%odi%nq2-1
                DO i = 1,vacuum%nmzxy
                   mapvac2 = mapvac2 + 1
                   j = j + 1
                   sout(j) =  AIMAG(den%vacxy(i,k,iv,3))
                END DO
             END DO
          END DO
          nvaccoeff2 = 2*vacuum%nmzxy*(oneD%odi%nq2-1)*vacuum%nvac + 2*vacuum%nmz*vacuum%nvac
          IF (mapvac2 .NE. nvaccoeff2) THEN
             WRITE (6,*)'The number of vaccum coefficients off the'
             WRITE (6,*)'off-diagonal part of the density matrix is'
             WRITE (6,*)'inconsitent:'
             WRITE (6,8000) mapvac2,nvaccoeff2
8000         FORMAT ('mapvac2= ',i12,'nvaccoeff2= ',i12)
             CALL juDFT_error("brysh1:# of vacuum coeff. inconsistent" ,calledby ="brysh1")
          ENDIF
       END IF
    ENDIF ! noco

    IF (atoms%n_u > 0 ) THEN     ! lda+U
       DO js = 1,input%jspins
          DO n = 1, atoms%n_u
             DO k = -3, 3
                DO i = -3, 3
                   j = j + 1 
                   sout(j) = REAL(den%mmpMat(i,k,n,js))
                   j = j + 1 
                   sout(j) = AIMAG(den%mmpMat(i,k,n,js))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF (input%film) THEN
       nvaccoeff = vacfac*vacuum%nmzxy*(oneD%odi%nq2-1)*vacuum%nvac + vacuum%nmz*vacuum%nvac
       IF (mapvac .NE. nvaccoeff) THEN
          WRITE(6,*)'The number of vaccum coefficients is'
          WRITE(6,*)'inconsitent:'
          WRITE (6,8010) mapvac,nvaccoeff
8010      FORMAT ('mapvac= ',i12,'nvaccoeff= ',i12)
          CALL juDFT_error("brysh1: # of vacuum coeff. inconsistent" ,calledby ="brysh1")
       ENDIF
    ENDIF

    mapmtd = atoms%ntype*(sphhar%nlhd+1)*atoms%jmtd
    IF (mapmt .GT. mapmtd) THEN
       WRITE(6,*)'The number of mt coefficients is larger than the'
       WRITE(6,*)'dimensions:'
       WRITE (6,8040) mapmt,mapmtd
8040   FORMAT ('mapmt= ',i12,' > mapmtd= ',i12)
       CALL juDFT_error("brysh1: mapmt > mapmtd (dimensions)",calledby ="brysh1")
    ENDIF

    nmap = j
    nall = (intfac*stars%ng3 + mapmt + mapvac + 49*2*atoms%n_u )*input%jspins
    IF (noco%l_noco) nall = nall + 2*stars%ng3 + mapvac2
    IF (nall.NE.nmap) THEN
       WRITE(6,*)'The total number of charge density coefficients is'
       WRITE(6,*)'inconsitent:'
       WRITE (6,8020) nall,nmap
8020   FORMAT ('nall= ',i12,'not equal nmap= ',i12)
       WRITE (6,'(a,i5,a,i5)') 'nall = ',nall,' nmap = ',nmap
       CALL juDFT_error ("brysh1: input # of charge density coeff. inconsistent" ,calledby ="brysh1")
    ENDIF
    IF (nmap.GT.SIZE(sout)) THEN 
       WRITE(6,*)'The total number of charge density coefficients is'
       WRITE(6,*)'larger than the dimensions:'
       WRITE (6,8030) nmap,SIZE(sout)
8030   FORMAT ('nmap= ',i12,' > size(sout)= ',i12)
       CALL juDFT_error("brysh1: nmap > mmap (dimensions)",calledby ="brysh1")
    ENDIF

  END SUBROUTINE brysh1
END MODULE m_brysh1
