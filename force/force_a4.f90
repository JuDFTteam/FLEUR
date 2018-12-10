MODULE m_force_a4
CONTAINS
  SUBROUTINE force_a4(atoms,sphhar,input,dimension,&
       &                    vr,&
       &                    force)
    !
    ! ************************************************************
    ! rhocore force contribution a la Rici et al.
    !
    ! ************************************************************
    !
    USE m_intgr, ONLY : intgr0,intgr3
    USE m_constants, ONLY : sfp_const
    USE m_differentiate,ONLY: difcub
    USE m_types
    USE m_cdn_io
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_dimension),INTENT(IN) :: dimension
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
    REAL,    INTENT (INOUT) :: force(3,atoms%ntype,input%jspins)
    !     ..
    !     .. Local Scalars ..
    REAL a4_1,a4_2,qcore,s13,s23,w,xi
    INTEGER i,ir,jsp,lh,lindex,mem,mindex,n,nd,na
    !     ..
    !     .. Local Arrays ..
    COMPLEX forc_a4(3),gv(3),ycomp1(3,-1:1)
    REAL rhoaux(atoms%jmtd),rhoc(atoms%jmtd,atoms%ntype,input%jspins)
    REAL tec(atoms%ntype,input%jspins),qintc(atoms%ntype,input%jspins)
    !     ..
    !     .. Data statements ..
    COMPLEX,PARAMETER:: czero=(0.000,0.000)
    !     ..
    !
    !     set ylm related components
    !
    s13 = SQRT(1.0/3.0)
    s23 = SQRT(2.0/3.0)
    ycomp1(1,0) = czero
    ycomp1(2,0) = czero
    ycomp1(3,0) = CMPLX(2.0*s13,0.0)
    ycomp1(1,-1) = CMPLX(s23,0.0)
    ycomp1(2,-1) = CMPLX(0.0,-s23)
    ycomp1(3,-1) = czero
    ycomp1(1,1) = CMPLX(-s23,0.0)
    ycomp1(2,1) = CMPLX(0.0,-s23)
    ycomp1(3,1) = czero
    !     --->    read in core density
    CALL readCoreDensity(input,atoms,dimension,rhoc,tec,qintc)

    DO jsp = 1,input%jspins
       na = 1
       DO n = 1,atoms%ntype
          IF (atoms%l_geo(n)) THEN
             DO i = 1,3
                forc_a4(i) = czero
             END DO
             !
             CALL intgr0(rhoc(:,n,jsp),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),qcore)
             !     write(16,1616) qcore
8000         FORMAT (' FORCE_A4: core charge=',1p,e16.8)
             !
             !
             nd = atoms%ntypsy(na)
             !
             !
             lh_loop: DO  lh = 1,sphhar%nlh(nd)
                lindex = sphhar%llh(lh,nd)
                IF (lindex.GT.1) EXIT lh_loop
                DO i = 1,3
                   gv(i) = czero
                END DO
                !
                !    sum over all m for particular symm. harmonic
                DO mem = 1,sphhar%nmem(lh,nd)
                   mindex = sphhar%mlh(mem,lh,nd)
                   DO i = 1,3
                      gv(i) = gv(i) + sphhar%clnu(mem,lh,nd)*ycomp1(i,mindex)
                   END DO
                END DO
                !
                !
                !     construct integrand rhocore*d/dr(v)*r**2
                !     note: rhocore is already multiplied by r**2 and srt(4.*pi)
                !     difcub performs analytic derivative of Lagrangian of 3rd order
                !
                xi = atoms%rmsh(1,n)
                rhoaux(1) = difcub(atoms%rmsh(1,n),vr(1,lh,n,jsp),xi)*rhoc(1,n,jsp)
                DO ir = 2,atoms%jri(n) - 2
                   xi = atoms%rmsh(ir,n)
                   rhoaux(ir) = difcub(atoms%rmsh(ir-1,n),&
                        &                                vr(ir-1,lh,n,jsp),xi) * rhoc(ir,n,jsp)
                END DO
                !
                ir = atoms%jri(n) - 1
                xi = atoms%rmsh(ir,n)
                rhoaux(ir) = difcub(atoms%rmsh(atoms%jri(n)-3,n),&
                     &                              vr(atoms%jri(n)-3,lh,n,jsp),xi)*rhoc(ir,n,jsp)
                !
                ir = atoms%jri(n)
                xi = atoms%rmsh(ir,n)
                rhoaux(ir) = difcub(atoms%rmsh(atoms%jri(n)-3,n),&
                     &                              vr(atoms%jri(n)-3,lh,n,jsp),xi)*rhoc(ir,n,jsp)
                CALL intgr3(rhoaux,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),w)
                a4_1 = 0.5*w/sfp_const
                !
                !     construct integrand rhocore*v*r
                !     note: rhocore is already multiplied by r**2 and srt(4.*pi)
                !
                DO ir = 1,atoms%jri(n)
                   rhoaux(ir) = rhoc(ir,n,jsp)/atoms%rmsh(ir,n)*vr(ir,lh,n,jsp)
                END DO
                !
                CALL intgr3(rhoaux,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),w)
                a4_2 = w/sfp_const
                !
                DO i = 1,3
                   forc_a4(i) = forc_a4(i) - (a4_1+a4_2)*gv(i)
                END DO
                !
                !  lh loop ends
             ENDDO lh_loop
             !
             !     sum to existing forces
             !
             DO i = 1,3
                force(i,n,jsp) = force(i,n,jsp) + REAL(forc_a4(i))
             END DO
             !
             !     write result
             !
             WRITE (6,FMT=8010) n
             WRITE (16,FMT=8010) n
             WRITE (6,FMT=8020) (forc_a4(i),i=1,3)
             WRITE (16,FMT=8020) (forc_a4(i),i=1,3)
8010         FORMAT (' FORCES: EQUATION A4 FOR ATOM TYPE',i4)
8020         FORMAT (' FX_A4=',2f10.6,' FY_A4=',2f10.6,' FZ_A4=',2f10.6)
             ! type loop ends
          ENDIF
          na = na + atoms%neq(n)
       ENDDO
       ! spin loop ends
    ENDDO
  END SUBROUTINE force_a4
END MODULE m_force_a4
