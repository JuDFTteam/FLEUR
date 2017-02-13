MODULE m_metric
  USE m_juDFT
  !*********************************************************
  !     multiplicates the vector s_in with the metric G
  !     output vector sout
  !********************************************************* 
CONTAINS
  SUBROUTINE metric(&
       &                  cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
       &                  mmap,nmaph,mapmt,mapvac2,s_in,sout,lpot) 
    USE m_metrz0
    USE m_convol
    USE m_types
    IMPLICIT NONE
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: mmap
    INTEGER, INTENT (IN) :: mapmt,mapvac2
    INTEGER, INTENT (IN) :: nmaph      
    LOGICAL, OPTIONAL,INTENT (IN) :: lpot !do we mix a potential??
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: s_in(mmap)
    REAL,    INTENT (OUT):: sout(mmap)
    !     ..
    !     .. Local Scalars ..
    INTEGER :: imap,ivac,iz,j,js,k2,l,n,iv2c,iv2,na,ioff
    REAL    :: dvol,dxn,dxn2,dxn4,volnstr2
    LOGICAL :: l_pot
    !     ..
    !     .. Local Arrays ..
    REAL,    ALLOCATABLE :: g(:),wght(:)
    COMPLEX, ALLOCATABLE :: ag3(:),fg3(:)
    !
    !     calculate the coefficients of the g-matrix 
    !     for the m.t. and vacuum region
    !
    IF (PRESENT(lpot)) THEN
       l_pot=lpot
    ELSE
       l_pot=.FALSE.
    ENDIF
    ALLOCATE (g(mmap),wght(vacuum%nmzd),ag3(stars%ng3),fg3(stars%ng3))
    g=0.0
    IF (sym%invs) THEN
       imap=stars%ng3
    ELSE
       imap=2*stars%ng3
    ENDIF
    iv2=1
    IF (.NOT.sym%invs2) iv2=2
    !
    ! metric for MT is r^2 dr/di di = r^3 dx ; divided by r^4 to 
    ! compensate use of n(r) * r^2 in array rho
    ! simpson integration used, weights for first and last point: 1
    ! weights forthe rest alternating: 2 or 4    
    !
    na = 1
    DO n = 1,atoms%ntype
       dxn = atoms%neq(n)*atoms%dx(n)/3.0e0
       dxn2 =2.0e0 *dxn
       dxn4 =4.0e0 *dxn
       !
       DO l = 0,sphhar%nlh(atoms%ntypsy(na))
          imap = imap + 1
          g(imap) = dxn/atoms%rmsh(1,n)
          IF (.NOT.l_pot) THEN
             DO j = 2,atoms%jri(n)-1,2
                imap = imap + 2
                g(imap-1) = dxn4/atoms%rmsh(j,n) 
                g(imap) = dxn2/atoms%rmsh(j+1,n) 
             ENDDO
             ! CHANGE JR 96/12/01
             ! take care when jri(n) is even
             imap=imap+1-MOD(atoms%jri(n),2)
             g(imap) = dxn/atoms%rmsh(atoms%jri(n),n)
          ELSE
             !
             !             for the potential multiply by r^4
             !
             DO j = 2,atoms%jri(n)-1,2
                imap = imap + 2
                g(imap-1) = dxn4*atoms%rmsh(j,n)**3 
                g(imap) = dxn2*atoms%rmsh(j+1,n)**3
             ENDDO
             imap=imap+1-MOD(atoms%jri(n),2)
             g(imap) = dxn*atoms%rmsh(atoms%jri(n),n)**3
          ENDIF

       ENDDO
       na = na + atoms%neq(n)
    enddo
    !
    ! vacuum contribution
    !
    IF (input%film) THEN
       dvol = cell%area*vacuum%delz
       !     nvac=1 if (zrfs.or.invs)
       IF ( vacuum%nvac.EQ.1 ) dvol = dvol + dvol
       IF (oneD%odi%d1) THEN
          dvol = cell%area*vacuum%delz
       END IF
       DO  ivac = 1,vacuum%nvac
          !                                      G||=0 components
          !
          !---> use 7-point simpson integration in accordance to intgz0.f
          !     calculate weights for integration
          !
          CALL metr_z0(vacuum%nmz,wght)
          DO iz = 1,vacuum%nmz
             imap = imap + 1
             IF (oneD%odi%d1) THEN
                g(imap) = wght(iz)*dvol*(cell%z1+(iz-1)*vacuum%delz)
             ELSE
                g(imap) = wght(iz)*dvol
             ENDIF
          ENDDO
          !                                      G||.ne.0 components
          !     calculate weights for integration
          CALL metr_z0(vacuum%nmzxy,wght)


          DO iv2c=1,iv2
             DO k2 = 1,oneD%odi%nq2-1
                IF (oneD%odi%d1) THEN
                   DO iz = 1,vacuum%nmzxy
                      imap = imap + 1
                      g(imap) = wght(iz)*oneD%odi%nst2(k2)*&
                           &                         dvol*(cell%z1+(iz-1)*vacuum%delz)
                   ENDDO
                ELSE
                   volnstr2= dvol*stars%nstr2(k2)
                   DO iz = 1,vacuum%nmzxy
                      imap = imap + 1
                      g(imap) = wght(iz)*volnstr2
                   ENDDO
                END IF
             ENDDO
          ENDDO
       enddo
    END IF
    !
    !     multiplicate the metric with the vector
    !
    DO  js = 1,input%jspins
       !    map s_in on a complex help array ag3
       IF (sym%invs) THEN
          DO imap = 1,stars%ng3
             ag3(imap) = CMPLX(s_in(imap+nmaph*(js-1)),0.0)
          ENDDO
       ELSE
          DO imap = 1,stars%ng3
             ag3(imap) = CMPLX(s_in(imap+nmaph*(js-1)),&
                  &                        s_in(imap+stars%ng3+nmaph*(js-1)))
          ENDDO
       ENDIF
       CALL convol(&
            &               stars,&
            &               fg3,&
            &               ag3)
       IF (sym%invs) THEN
          DO imap = 1,stars%ng3
             sout(imap+nmaph*(js-1)) = cell%omtil*REAL(fg3(imap))
          ENDDO
          DO imap = stars%ng3+1,nmaph
             sout(imap+nmaph*(js-1)) = g(imap)*s_in(imap+nmaph*(js-1))
          ENDDO
       ELSE
          DO imap = 1,stars%ng3
             sout(imap+nmaph*(js-1)) = cell%omtil*REAL(fg3(imap))
             sout(imap+stars%ng3+nmaph*(js-1)) = cell%omtil*AIMAG(fg3(imap))
          ENDDO
          DO imap = 2*stars%ng3+1,nmaph
             sout(imap+nmaph*(js-1)) = g(imap)*s_in(imap+nmaph*(js-1))
          ENDDO
       ENDIF
    enddo

    IF (noco%l_noco) THEN
       DO imap = 1,stars%ng3
          ag3(imap) = CMPLX(s_in(2*nmaph + imap),&
               &                        s_in(2*nmaph + stars%ng3 + imap))
       ENDDO
       CALL convol(&
            &        stars,&
            &        fg3,&
            &        ag3)
       DO imap = 1,stars%ng3
          sout(2*nmaph + imap)       = cell%omtil*REAL(fg3(imap))
          sout(2*nmaph + stars%ng3 + imap) = cell%omtil*AIMAG(fg3(imap))
       ENDDO
       IF (input%film) THEN
          !--->    js runs over the real and imaginary part of the vacuum density
          !--->    coefficients (not the spin).
          IF (sym%invs) THEN
             ioff = stars%ng3 + mapmt
          ELSE
             ioff = 2*stars%ng3 + mapmt
          ENDIF
          DO js = 1,2
             DO imap = 1,mapvac2/2
                sout(2*nmaph + 2*stars%ng3 + mapvac2/2*(js-1) + imap) = &
                     &               g(ioff + imap)&
                     &               *s_in(2*nmaph + 2*stars%ng3 + mapvac2/2*(js-1) + imap)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    IF ( atoms%n_u > 0 )  THEN
       j = input%jspins*nmaph
       IF (noco%l_noco) THEN
          j = j +  2*stars%ng3 
          IF (input%film) THEN
             j = j + mapvac2
          ENDIF
       ENDIF
       DO imap = j+1, j+49*2*input%jspins*atoms%n_u
          sout(imap) = s_in(imap)
       ENDDO
    ENDIF
    DEALLOCATE (g,wght,ag3,fg3)

  END SUBROUTINE metric
END MODULE m_metric
