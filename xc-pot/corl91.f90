    MODULE m_corl91
!.....-----------------------------------------------------------------
!     uniform-gas correlation of perdew and wang 1991
!.....-----------------------------------------------------------------
    CONTAINS
    SUBROUTINE corl91( &
    rs,zta, &
    ec,vcup,vcdn,ecrs,eczta,alfc)
!.....-----------------------------------------------------------------
!     input: seitz radius (rs), relative spin polarization (zta)
!     output: correlation energy per electron (ec),
!             up- and down-spin potentials (vcup,vcdn),
!             derivatives of ec wrt rs (ecrs) &zta (eczta).
!     output: correlation contribution (alfc) to the spin stiffness
!.....-----------------------------------------------------------------
    IMPLICIT NONE

    REAL, INTENT (IN)  :: rs,zta
    REAL, INTENT (OUT) :: alfc,ec,ecrs,eczta,vcdn,vcup

    REAL :: alfm,alfrsm,comm,ep,eprs,eu,eurs,f,fz,z4
    REAL :: fzz,gam,thrd,thrd4
!.....-----------------------------------------------------------------
!     ..
    DATA gam,fzz/0.5198421,1.709921/
    DATA thrd,thrd4/0.333333333333e0,1.333333333333e0/
!.....-----------------------------------------------------------------
    f = ((1.0+zta)**thrd4+ (1.0-zta)**thrd4-2.e0)/gam

    CALL gcor91(0.0310907,0.21370,7.5957,3.5876,1.6382, &
    &           0.49294,1.00,rs,eu,eurs)

    CALL gcor91(0.01554535,0.20548,14.1189,6.1977,3.3662, &
    &           0.62517,1.00,rs,ep,eprs)

    CALL gcor91(0.0168869,0.11125,10.357,3.6231,0.88026, &
    &           0.49671,1.00,rs,alfm,alfrsm)

!  alfm is minus the spin stiffness alfc
    alfc = -alfm
    z4 = zta**4
    ec = eu* (1.0-f*z4) + ep*f*z4 - alfm*f* (1.0-z4)/fzz
!  energy done. now the potential:
    ecrs = eurs* (1.0-f*z4) + eprs*f*z4 - alfrsm*f* (1.0-z4)/fzz
    fz = thrd4* ((1.0+zta)**thrd- (1.0-zta)**thrd)/gam
    eczta = 4.e0* (zta**3)*f* (ep-eu+alfm/fzz) + &
    fz* (z4*ep-z4*eu- (1.0-z4)*alfm/fzz)
    comm = ec - rs*ecrs/3.e0 - zta*eczta
    vcup = comm + eczta
    vcdn = comm - eczta

    END SUBROUTINE corl91
!.....-----------------------------------------------------------------
!     called by corl91
!.....-----------------------------------------------------------------
    SUBROUTINE gcor91( &
    a,a1,b1,b2,b3,b4,p,rs, &
    gg,ggrs)

    IMPLICIT NONE
    REAL, INTENT (IN)  :: a,a1,b1,b2,b3,b4,p,rs
    REAL, INTENT (OUT) :: gg,ggrs

    REAL :: p1,q0,q1,q2,q3,rs12,rs32,rsp
!.....-----------------------------------------------------------------
!     ..
    p1 = p + 1.0
    q0 = -2.e0*a* (1.0+a1*rs)
    rs12 = sqrt(rs)
    rs32 = rs12**3
    rsp = rs**p
    q1 = 2.e0*a* (b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
    q2 = log(1.0+1.0/q1)
    gg = q0*q2
    q3 = a* (b1/rs12+2.e0*b2+3.e0*b3*rs12+2.e0*b4*p1*rsp)
    ggrs = -2.e0*a*a1*q2 - q0*q3/ (q1**2+q1)

    END SUBROUTINE gcor91

    END MODULE m_corl91
