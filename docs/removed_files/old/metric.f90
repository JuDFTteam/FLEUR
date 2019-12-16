MODULE m_metric
  USE m_juDFT
  !*********************************************************
  !     multiplicates the vector s_in with the metric G
  !     output vector sout
  !********************************************************* 
CONTAINS
  SUBROUTINE metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
                    mmap,nmaph,mapmt,mapvac2,s_in,sout,lpot) 
    USE m_metrz0
    USE m_convol
    USE m_types

    IMPLICIT NONE

    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_vacuum),INTENT(IN) :: vacuum
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_stars),INTENT(IN)  :: stars
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_sphhar),INTENT(IN) :: sphhar
    TYPE(t_atoms),INTENT(IN)  :: atoms

    ! Scalar Arguments
    INTEGER, INTENT (IN)          :: mmap
    INTEGER, INTENT (IN)          :: mapmt,mapvac2
    INTEGER, INTENT (IN)          :: nmaph      
    LOGICAL, OPTIONAL,INTENT (IN) :: lpot !do we mix a potential??

    ! Array Arguments
    REAL,    INTENT (IN)          :: s_in(mmap)
    REAL,    INTENT (OUT)         :: sout(mmap)

    ! Local Scalars
    INTEGER              :: imap,ivac,iz,j,js,k2,l,n,iv2c,iv2,na,ioff,i
    REAL                 :: dvol,dxn,dxn2,dxn4,volnstr2
    LOGICAL              :: l_pot

    ! Local Arrays
    REAL,    ALLOCATABLE :: g(:),wght(:)
    COMPLEX, ALLOCATABLE :: ag3(:),fg3(:)

    ! calculate the coefficients of the g-matrix
    ! for the m.t. and vacuum region

    l_pot = .FALSE.
    IF (PRESENT(lpot)) l_pot = lpot ! for potential mixing

    ALLOCATE (g(mmap),wght(vacuum%nmzd),ag3(stars%ng3),fg3(stars%ng3))
    g = 0.0

    imap = 2 * stars%ng3 ! complex values without invs
    IF (sym%invs) imap = stars%ng3 ! only real values with invs
    iv2 = 2
    IF (sym%invs2) iv2 = 1

    ! metric for MT is r^2 dr/di di = r^3 dx ; divided by r^4 to 
    ! compensate use of n(r) * r^2 in array rho
    ! simpson integration used, weights for first and last point: 1
    ! weights forthe rest alternating: 2 or 4    
    na = 1
    DO n = 1, atoms%ntype
       dxn = atoms%neq(n) * atoms%dx(n) / 3.0e0
       dxn2 =2.0e0 * dxn
       dxn4 =4.0e0 * dxn

       DO l = 0, sphhar%nlh(sym%ntypsy(na))
          imap = imap + 1
          g(imap) = dxn / atoms%rmsh(1,n)
          IF (.NOT.l_pot) THEN
             DO j = 2, atoms%jri(n) - 1, 2
                imap = imap + 2
                g(imap-1) = dxn4 / atoms%rmsh(j,n) 
                g(imap) = dxn2 / atoms%rmsh(j+1,n) 
             END DO
             ! CHANGE JR 96/12/01
             ! take care when jri(n) is even
             imap = imap + 1 - MOD(atoms%jri(n),2)
             g(imap) = dxn / atoms%rmsh(atoms%jri(n),n)
          ELSE
             ! for the potential multiply by r^4
             DO j = 2, atoms%jri(n) - 1, 2
                imap = imap + 2
                g(imap-1) = dxn4 * atoms%rmsh(j,n)**3 
                g(imap) = dxn2 * atoms%rmsh(j+1,n)**3
             END DO
             imap = imap + 1 - MOD(atoms%jri(n),2)
             g(imap) = dxn * atoms%rmsh(atoms%jri(n),n)**3
          END IF
       END DO
       na = na + atoms%neq(n)
    END DO

    ! vacuum contribution
    IF (input%film) THEN
       dvol = cell%area*vacuum%delz
       ! nvac=1 if (zrfs.or.invs)
       IF (vacuum%nvac.EQ.1) dvol = dvol + dvol
       IF (oneD%odi%d1) dvol = cell%area*vacuum%delz
       DO ivac = 1, vacuum%nvac
          ! G||=0 components
          !
          ! use 7-point simpson integration in accordance to intgz0.f
          ! calculate weights for integration
          CALL metr_z0(vacuum%nmz,wght)
          DO iz = 1, vacuum%nmz
             imap = imap + 1
             IF (oneD%odi%d1) THEN
                g(imap) = wght(iz) * dvol * (cell%z1+(iz-1)*vacuum%delz)
             ELSE
                g(imap) = wght(iz) * dvol
             END IF
          END DO
          ! G||.ne.0 components
          !
          ! calculate weights for integration
          CALL metr_z0(vacuum%nmzxy,wght)
          DO iv2c = 1, iv2
             DO k2 = 1, oneD%odi%nq2 - 1
                IF (oneD%odi%d1) THEN
                   DO iz = 1,vacuum%nmzxy
                      imap = imap + 1
                      g(imap) = wght(iz) * oneD%odi%nst2(k2) * dvol * (cell%z1+(iz-1)*vacuum%delz)
                   END DO
                ELSE
                   volnstr2 = dvol * stars%nstr2(k2)
                   DO iz = 1, vacuum%nmzxy
                      imap = imap + 1
                      g(imap) = wght(iz) * volnstr2
                   END DO
                END IF
             END DO
          END DO
       END DO
    END IF

    ! Apply metric to interstitial region (metric here = step function)
    ! + multiply metric g with s_in for MT and vacuum contributions (store in sout)
    DO js = 1, input%jspins
       ! map s_in on a complex help array ag3
       IF (sym%invs) THEN
          DO imap = 1, stars%ng3
             ag3(imap) = CMPLX(s_in(imap+nmaph*(js-1)),0.0)
          END DO
       ELSE
          DO imap = 1, stars%ng3
             ag3(imap) = CMPLX(s_in(imap+nmaph*(js-1)),s_in(imap+stars%ng3+nmaph*(js-1)))
          END DO
       ENDIF

       CALL convol(stars,fg3,ag3,stars%ufft)

       IF (sym%invs) THEN
          ! interstitial
          DO imap = 1, stars%ng3
             sout(imap+nmaph*(js-1)) = cell%omtil * REAL(fg3(imap))
          END DO
          ! MT + vacuum
          DO imap = stars%ng3+1, nmaph
             sout(imap+nmaph*(js-1)) = g(imap) * s_in(imap+nmaph*(js-1))
          END DO
       ELSE
          ! interstitial
          DO imap = 1, stars%ng3
             sout(imap+nmaph*(js-1))           = cell%omtil * REAL(fg3(imap))
             sout(imap+stars%ng3+nmaph*(js-1)) = cell%omtil * AIMAG(fg3(imap))
          END DO
          ! MT + vacuum
          DO imap = 2 * stars%ng3 + 1, nmaph
             sout(imap+nmaph*(js-1)) = g(imap) * s_in(imap+nmaph*(js-1))
          END DO
       END IF
    END DO

    IF (noco%l_noco) THEN
       DO imap = 1, stars%ng3
          ag3(imap) = CMPLX(s_in(2*nmaph + imap),s_in(2*nmaph + stars%ng3 + imap))
       END DO

       CALL convol(stars,fg3,ag3,stars%ufft)

       DO imap = 1, stars%ng3
          sout(2*nmaph + imap)             = cell%omtil * REAL(fg3(imap))
          sout(2*nmaph + stars%ng3 + imap) = cell%omtil * AIMAG(fg3(imap))
       END DO
       IF (input%film) THEN
          ! js runs over the real and imaginary part of the vacuum density
          ! coefficients (not the spin).
          ioff = 2*stars%ng3 + mapmt ! (for complex values)
          IF (sym%invs) ioff = stars%ng3 + mapmt ! (real values if sym%invs)

          DO js = 1, 2
             DO imap = 1, mapvac2 / 2
                sout(2*nmaph + 2*stars%ng3 + mapvac2/2*(js-1) + imap) = &
                   g(ioff + imap) * s_in(2*nmaph + 2*stars%ng3 + mapvac2/2*(js-1) + imap)
             END DO
          END DO
       END IF
       !mt offdiagonal part
       imap=2*nmaph + 2*stars%ng3 + mapvac2/2*(js-1) + imap-1
       ioff=MERGE(stars%ng3,2*stars%ng2,sym%invs)+1
       j=0
       IF (noco%l_mtnocopot) THEN
          na = 1
          DO n = 1,atoms%ntype
             DO l = 0,sphhar%nlh(sym%ntypsy(na))
                DO i = 1,atoms%jri(n)
                   j = j + 1
                   sout(imap+j) = g(ioff+(j-1)/2)*s_in(imap+j)
                   j = j + 1
                   sout(imap+j) = g(ioff+(j-1)/2)*s_in(imap+j)
                END DO
             END DO
             na = na + atoms%neq(n)
          END DO
       END IF
    END IF

    ! density matrix
    IF (atoms%n_u > 0) THEN
       j = input%jspins * nmaph
       IF (noco%l_noco) THEN
          j = j + 2 * stars%ng3
          IF (input%film) THEN
             j = j + mapvac2
          END IF
       END IF
       DO imap = j+1, j+49*2*input%jspins*atoms%n_u
          sout(imap) = s_in(imap)
       END DO
    END IF

    DEALLOCATE (g,wght,ag3,fg3)

  END SUBROUTINE metric
END MODULE m_metric
