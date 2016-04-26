MODULE m_brysh2
  USE m_juDFT
  !******************************************************
  !     maps the density back from one single vector into the
  !     proper component of interstitial, m.t. and vacuum density
  !******************************************************
CONTAINS 
  SUBROUTINE brysh2(&
       &                  input,stars,atoms,sphhar,&
       &                  noco,vacuum,&
       &                  sym,s_in,&
       &                  n_mmp,oneD,qpw,rho,rht,rhtxy,cdom,cdomvz,cdomvxy)
    USE m_types
    IMPLICIT NONE

    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    REAL,    INTENT (IN) :: s_in(:)
    REAL,    INTENT (OUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins)
    REAL,    INTENT (OUT) :: rht(vacuum%nmz,2,input%jspins)
    COMPLEX, INTENT (OUT) :: qpw(stars%n3d,input%jspins),cdom(stars%n3d),cdomvz(vacuum%nmz,2)
    COMPLEX, INTENT (OUT) :: rhtxy(vacuum%nmzxy,oneD%odi%n2d-1,2,input%jspins)
    COMPLEX, INTENT (OUT) :: cdomvxy(vacuum%nmzxy,oneD%odi%n2d-1,2)
    COMPLEX, INTENT (OUT) :: n_mmp(-3:3,-3:3,atoms%n_u,input%jspins)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,iv,j,js,k,l,n,na
    !
    j=0
    DO  js = 1,input%jspins
       IF (sym%invs) THEN
          DO i = 1,stars%ng3
             j = j + 1
             qpw(i,js) = CMPLX(s_in(j),0.0)
          END DO
       ELSE
          DO i = 1,stars%ng3
             j = j + 1
             qpw(i,js) = CMPLX(s_in(j),s_in(j+stars%ng3))
          END DO
          j = j + stars%ng3
       ENDIF
       na = 1
       DO n = 1,atoms%ntype
          DO l = 0,sphhar%nlh(atoms%ntypsy(na))
             DO i = 1,atoms%jri(n)
                j = j + 1
                rho(i,l,n,js) = s_in(j)
             END DO
          END DO
          na = na + atoms%neq(n)
       END DO
       IF (input%film) THEN
          DO iv = 1,vacuum%nvac
             DO k = 1,vacuum%nmz
                j = j + 1
                rht(k,iv,js) = s_in(j)
             END DO
             DO k = 1,oneD%odi%nq2-1
                DO i = 1,vacuum%nmzxy
                   j = j + 1
                   rhtxy(i,k,iv,js) = CMPLX(s_in(j),0.0)
                END DO
             END DO
             IF (.NOT.sym%invs2) THEN
                DO k = 1,oneD%odi%nq2-1
                   DO i = 1,vacuum%nmzxy
                      j = j + 1
                      rhtxy(i,k,iv,js) = rhtxy(i,k,iv,js) +CMPLX(0.0,s_in(j))
                   END DO
                END DO
             END IF
          END DO
       END IF
    enddo

    IF (noco%l_noco) THEN
       !--->    off-diagonal part of the density matrix
       DO i = 1,stars%ng3
          j = j + 1
          cdom(i) = CMPLX(s_in(j),0.0)
       END DO
       DO i = 1,stars%ng3
          j = j + 1
          cdom(i) = cdom(i) + CMPLX(0.0,s_in(j))
       END DO
       IF (input%film) THEN
          DO iv = 1,vacuum%nvac
             DO k = 1,vacuum%nmz
                j = j + 1
                cdomvz(k,iv) = CMPLX(s_in(j),0.0)
             END DO
             DO k = 1,oneD%odi%nq2-1
                DO i = 1,vacuum%nmzxy
                   j = j + 1
                   cdomvxy(i,k,iv) = CMPLX(s_in(j),0.0)
                END DO
             END DO
          END DO
          DO iv = 1,vacuum%nvac
             DO k = 1,vacuum%nmz
                j = j + 1
                cdomvz(k,iv) = cdomvz(k,iv) + CMPLX(0.0,s_in(j))
             END DO
             DO k = 1,oneD%odi%nq2-1
                DO i = 1,vacuum%nmzxy
                   j = j + 1
                   cdomvxy(i,k,iv) = cdomvxy(i,k,iv)+ CMPLX(0.0,s_in(j))
                END DO
             END DO
          END DO
       END IF
    ENDIF

    IF ( atoms%n_u > 0 ) THEN
       DO js = 1,input%jspins
          DO n = 1, atoms%n_u
             DO k = -3, 3
                DO i = -3, 3
                   j = j + 1
                   n_mmp(i,k,n,js) = CMPLX(s_in(j),s_in(j+1))
                   j = j + 1
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE brysh2
   END MODULE m_brysh2
