!--------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

 MODULE m_differ
    use m_juDFT
!
!differ   subroutine for the differential equations
!-----x----x----x----x----x----x----x----x----x----x----x----x----x----
!-----x this version modified to accept energies up to emax=2     x----
!-----x                              d.d.koelling   7/7/77      x----x
!-----x----x----x----x----x----x----x----x----x----x----x----x----x----
    CONTAINS
    SUBROUTINE differ(&
                     fn,fl,fj,c,z,h,rnot,rn,d,msh,vr,&
                     e,&
                     a,b,ierr)
    USE m_constants
    USE m_inwint
    USE m_outint

    IMPLICIT NONE

!     .. Scalar Arguments ..
    INTEGER,INTENT(IN)  :: msh
    INTEGER,INTENT(OUT) :: ierr
    REAL,   INTENT(IN)  :: fj,fl,fn,c,z,h,rnot,rn,d
    REAL,   INTENT(IN)  :: vr(msh)
    REAL,   INTENT(OUT) :: a(msh),b(msh)
    REAL,   INTENT(INOUT) :: e

!     .. Local Scalars ..
    REAL cis,cs,de,del,dg,eabsv,emax,emin,fkap,g,qcoef,qqqq,r,ra,rb,rg,rj,s,w,wmin
    INTEGER k,ki,kj,n,nodes,nqnt,ntimes
    LOGICAL dbl
    CHARACTER(LEN=150) hintString

!     .. Local Arrays ..
    REAL a0(5),b0(5)

    !init 
    a = 0.0
    b = 0.0
    ierr = 0
    nqnt = fn - fl - 0.99e0
    n = msh
    del = 1.e-7
    emax = 2
    emin = -z*z/fn**2 - 10.0
    s      = 2.0e0* (fj-fl)
    cs     = c*s
    cis    = 1.0e0/cs
    fkap   = fj + 0.5e0
    g      = sqrt(fkap**2- (cis*z)**2)
    
    try_loop:DO ntimes=1,200
        !bisection
        IF ((e.GE.emax) .OR. (e.LE.emin)) e = 0.5e0* (emax+emin)
        !Do outward integration ....
        CALL outint(msh,e,fkap,cs,cis,s,vr,z,rn,rnot,h,d,a0,b0,a,b,ki,nodes)
        IF (nqnt-nodes.GT.0) THEN
            !**** too few nodes
            IF (e.GT.emin) emin = e
        ELSEIF (nqnt-nodes.LT.0) THEN
            !**** too many nodes
            IF (e.LT.emax) emax = e
        ELSEIF (nqnt-nodes.EQ.0) THEN
            !**** correct number of nodes
            ra = a(ki)
            rb = b(ki)
            CALL inwint(e,fl,ki,fkap,cs,cis,s,z,h,d,rn,rnot,msh,vr,a,b,kj)
            ra = ra/a(ki)
            rb = rb/b(ki)
            a(ki:kj)=a(ki:kj)*ra
            b(ki:kj)=b(ki:kj)*ra
            dg = exp(h*g)
            rg = rnot**g
            DO k = 1,kj
                a(k) = a(k)*rg
                b(k) = rg*b(k)
                rg   = rg*dg
            ENDDO
            eabsv = abs(e)
            qcoef = sqrt(eabsv+eabsv)
            rj = rnot*d** (kj-1)
            r = rj
            rg = a(kj)
            dg = b(kj)
            qqqq = min(abs(rg),abs(dg))
            IF (qqqq.GE.1.0e-25) THEN
                wmin = (1.e-35)/qqqq
                DO k=kj+1,n
                    r = r*d
                    w =  exp(qcoef* (rj-r))
                    IF (w.LT.wmin) exit
                    a(k) = w*rg
                    b(k) = w*dg
                enddo
            ENDIF
            a(k+1:n) = 0.0e0
            b(k+1:n) = 0.0e0
            r = rn
            w = r* (a(n)**2+b(n)**2)
            
            r = r + r
            rj = 1.0e0/d
            dbl = .false.
            DO k=n-1,1,-1
                r = r*rj
                rg = r* (a(k)**2+b(k)**2)
                w = w + rg
                dbl = .NOT. dbl
                IF (dbl) w = w + rg
            ENDDO
            
            w = h/3.0* (w+rnot* (a(1)**2+b(1)**2))
            de = cs*a(ki)*b(ki)* (ra-rb)/ (ra*w)   
            IF (de.GT.0.0e0) emin = e
            IF (de.LT.0.0e0) emax = e
            e = e + de
            IF (abs(de).LT.del) then
              
                exit try_loop 
            endif      
        endif
        !Check if interval is empty
        IF ((emax-emin).LE.del) THEN
            WRITE (oUnit,FMT=8010) nodes,fn,fl,fj,emin,e,emax
            WRITE (oUnit,FMT=8030) vr
            WRITE (oUnit,FMT=8030) a
            WRITE (hintString,'(a,i0,a,i0,a,i0,a)') "The n=",NINT(fn)," l=",NINT(fl), " state of an atom with Z=", NINT(z),&
                " seems to be not in a reasonable energy range."
            CALL juDFT_error("differ 1: problems with solving dirac equation",calledby ="differ", hint=TRIM(hintString))
        ENDIF    
    
    ENDDO try_loop

    !
    ! for a consistent definition of the small component in
    ! the Dirac and scalar-relativistic approximation (SRA)
    ! we have to multiply the small component with -1;
    ! then the small component of the Dirac equation
    ! also fulfills the SRA for s-states
    !                                               M.B. (July, 2012)
    a=a/sqrt(w)
    b=-b/sqrt(w)

    if (ntimes.gt.200) ierr=1 !we did not converge

    IF (ntimes.GT.200.and.(nqnt.ne.nodes)) THEN
        WRITE (oUnit,FMT=8000)
        WRITE (oUnit,FMT=8030) fn,fl,fj,emin,emax
        WRITE (oUnit,FMT=8030) e,de
        WRITE (oUnit,FMT=8030) ra,rb,w,a(ki),b(ki)
        WRITE (oUnit,FMT=8000)
        WRITE (oUnit,FMT=8030) vr
        WRITE (oUnit,FMT=8030) a
        CALL juDFT_error("Differ did not converge")
    8000 FORMAT (/,/,/,/,10x,' too many tries required')
    8010 FORMAT (/,/,/,/,10x,' too many nodes.',i5,3f4.1,3e15.7)
    8020 FORMAT (/,/,/,/,10x,' too few nodes. ',i5,3f4.1,3e15.7)
    8030 FORMAT (10x,5e14.4)
    ENDIF
    END SUBROUTINE differ
    END MODULE m_differ