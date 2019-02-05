MODULE m_force_a3
CONTAINS
  SUBROUTINE force_a3(atoms,sphhar,&
       &                    input,&
       &                    rho,vr,&
       &                    force)
    ! ************************************************************
    ! Hellman-Feynman force contribution a la Rici et al.
    ! ************************************************************
    !
    USE m_intgr, ONLY : intgr3
    USE m_constants, ONLY : pi_const
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     .. Array Arguments ..
    REAL,    INTENT (IN) ::  vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
    REAL,    INTENT (IN) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
    REAL,    INTENT (INOUT) :: force(3,atoms%ntype,input%jspins)
    !     ..
    !     .. Local Scalars ..
    COMPLEX,PARAMETER:: czero=(0.0,0.0)
    REAL a3_1,a3_2,s34,s38,w
    INTEGER i,ir,jsp,lh,lindex,mem,mindex,n,nd,na
    !     ..
    !     .. Local Arrays ..
    COMPLEX forc_a3(3),grd1(3,-1:1),gv(3)
    REAL rhoaux(atoms%jmtd)
    !     ..
    !     set components of gradient in terms of Ylms
    !
    s34 = SQRT(3.0/ (4.0*pi_const))
    s38 = SQRT(3.0/ (8.0*pi_const))
    grd1(1,0) = czero
    grd1(2,0) = czero
    grd1(3,0) = CMPLX(s34,0.0)
    grd1(1,-1) = CMPLX(s38,0.0)
    grd1(2,-1) = CMPLX(0.0,-s38)
    grd1(3,-1) = czero
    grd1(1,1) = CMPLX(-s38,0.0)
    grd1(2,1) = CMPLX(0.0,-s38)
    grd1(3,1) = czero
    !
    WRITE  (6,*)
    spin: DO jsp = 1,input%jspins
       na = 1
       DO n = 1,atoms%ntype
          IF (atoms%l_geo(n)) THEN

             DO i = 1,3
                forc_a3(i) = czero
             END DO
             !
             nd = atoms%ntypsy(na)
             !
             lat_har: DO lh = 1,sphhar%nlh(nd)
                lindex = sphhar%llh(lh,nd)
                IF (lindex.GT.1) EXIT
                DO i = 1,3
                   gv(i) = czero
                END DO
                !    sum over all m for particular symm. harmonic
                DO mem = 1,sphhar%nmem(lh,nd)
                   mindex = sphhar%mlh(mem,lh,nd)
                   DO i = 1,3
                      gv(i) = gv(i) + sphhar%clnu(mem,lh,nd)*grd1(i,mindex)
                   END DO
                END DO
                DO ir = 1,atoms%jri(n)
                   rhoaux(ir) = rho(ir,lh,n,jsp)/ (atoms%rmsh(ir,n)**2)*&
                        &                         (1.0- (atoms%rmsh(ir,n)/atoms%rmt(n))**3)
                END DO
                CALL intgr3(rhoaux,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),w)
                a3_1 = 4.0*pi_const/3.0*w
                !+Gu
                a3_2 = vr(atoms%jri(n),lh,n,jsp)/(input%jspins*atoms%rmt(n))
                !-Gu
                DO i = 1,3
                   forc_a3(i) = forc_a3(i) + (a3_1+a3_2)*gv(i)*atoms%zatom(n)
                END DO
             ENDDO lat_har

             DO i = 1,3
                force(i,n,jsp) = force(i,n,jsp) + REAL(forc_a3(i))
             END DO
             !
             !     write result
             WRITE (6,FMT=8010) n
             WRITE (6,FMT=8020) (forc_a3(i),i=1,3)
8010         FORMAT (' FORCES: EQUATION A3 FOR ATOM TYPE',i4)
8020         FORMAT (' FX_A3=',2f10.6,' FY_A3=',2f10.6,' FZ_A3=',2f10.6)

          ENDIF               ! atoms%l_geo(n) = .true.
          na = na + atoms%neq(n)
       ENDDO

    ENDDO spin

  END SUBROUTINE force_a3
END MODULE m_force_a3
