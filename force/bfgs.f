      MODULE m_bfgs
      use m_juDFT
      CONTAINS
      SUBROUTINE bfgs(
     >                ntype,istep,istep0,force,
     >                zat,xa,thetad,epsdisp,epsforce,tote,
     X                xold,y,h,tau0,
     <                lconv)
!*******************************************************************
!     relaxes the forces using the BFGS quasi-Newton method.
!     input:
!     istep  = atomic step in this run
!     istep0 = number of atomic steps in previous runs used
!              in updating hessian
!     
!     output:
!        lconv   = logical true if forces are converged to tolerance
!                  given in epsforce
!     
!     the positions and forces from this step are added 
!     to file force.dat
!*******************************************************************
      
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)    :: ntype,istep,istep0
      REAL,    INTENT (IN)    :: epsdisp,epsforce,tote
      REAL,    INTENT (INOUT) :: thetad,xa
      LOGICAL, INTENT (OUT)   :: lconv
!     ..
!     .. Array Arguments ..
      REAL,    INTENT (IN)    :: force(3,ntype),zat(ntype)
      REAL,    INTENT (INOUT) :: tau0(3,ntype)
      REAL,    INTENT (INOUT) :: y(3*ntype),xold(3*ntype)
      REAL,    INTENT (INOUT) :: h(3*ntype,3*ntype)
!     ..
!     .. Local Scalars ..
      INTEGER :: i,j,n,nn,ist,na
      REAL    :: py,yy,alpha,s,zatm,fmax,gamma,d2,dispmax
!     ..
!     .. Local  Arrays ..
      REAL :: f(3*ntype),p(3*ntype),v(3*ntype),xnew(3*ntype)
!     
      n=3*ntype
      ist = istep+istep0
      
!---  >       get positions and forces in correct form and output
      DO na = 1,ntype
         nn = 3*(na-1)
         xnew(nn+1:nn+3) = tau0(:,na)
         f(nn+1:nn+3)    = force(:,na)
      ENDDO
      !Write new entry into forces.dat
      OPEN(43,file ='forces.dat',status ='unknown',form='formatted'
     $     ,position ='append')
      WRITE (43,'(a,f20.10)') "energy =",tote
      WRITE (43,'(3f16.9,3f14.9)') 
     >     ((xnew(i+3*j),i = 1,3),(f(i+3*j),i = 1,3),j = 0,ntype-1)
      CLOSE(43)
!---  >       get maximum force
      fmax = 0.0
      DO na = 1,ntype
         nn = 3*(na-1)
         fmax = MAX( fmax, (f(nn+1)**2+f(nn+2)**2+f(nn+3)**2) )
      ENDDO
      fmax = SQRT(fmax)
      WRITE (6,1000) istep,fmax
      IF ( fmax<epsforce) THEN
         lconv = .TRUE.
         RETURN
      ELSE
         lconv = .FALSE.
      ENDIF
 1000 FORMAT (1x/,' atomic step',i4,': maximum force =',
     +     1p,e14.6,' hartrees/a.u.')

!------------------------------------------------------------
! ===>       if first step, go along gradient
      
      IF (ist==1) THEN
!---  >       choose a reasonable first guess for scaling, but
!---  >       limit displacement to a maximum of 0.25 a.u.
!---  >       (may need to be changed for different systems)
!---  >       this choice is based on a Debye temperature of 330K;
!---  >       modify as needed (change thetad in param.8)
         zatm = 0.0
         DO i = 1,ntype
            zatm = MAX(zatm,zat(i))
         ENDDO
         IF (ABS(xa)<1.0e-10) THEN
            WRITE (6,*) 'WARNING, xa = 0.0 set to 2.0'
            xa = 2.0
         ENDIF
         IF (ABS(thetad)<1.0e-10) THEN
            WRITE (6,*) 'WARNING, thetad = 0.0 set to 330.0'
            thetad = 330.0
         ENDIF
         alpha = (250.0/(zatm*xa))*((330./thetad)**2)
         IF ( alpha*fmax*xa > 0.15 ) alpha = 0.25/(fmax*xa)
         p(:) = alpha*f(:)
!---  >       set h to identity matrix
         h = 0.0
         DO j = 1,n
            h(j,j) = 1.0
         ENDDO
         
      ELSE
!------------------------------------------------------------
! ===>              k-th step
!     
!---  >       determine p
         p(:) = xnew(:)-xold(:)
!---  >       update the change in gradients
         y(:) = y(:)-f(:)
!---  >       get necessary inner products and H|y>
         py = dot_PRODUCT(p,y)
         yy = 0.0
         DO i = 1,n
            s = 0.0
            DO j = 1,n
               s = s+y(j)*h(j,i)
            ENDDO
            v(i) = s
            yy   = yy+y(i)*s
         ENDDO
!---  >       check that update will leave h positive definite;
!---  >       if not, then stop
         IF (py<=0.0) THEN
            WRITE (6,*) '  bfgs: <p|y> < 0'
            WRITE (6,*) '  check convergence of forces'
             CALL juDFT_error("bfgs: <p|y><0",calledby="bfgs")
         ELSE
!---  >          update h
            IF (ist==2) THEN
               gamma = py/yy
            ELSE
               gamma = 1.0
            ENDIF
            DO j = 1,n
               DO i = 1,n
                  h(i,j) = (h(i,j) - (v(i)*p(j)+p(i)*v(j))/py)*gamma
     +                 + (1.+gamma*yy/py)*p(i)*p(j)/py
               ENDDO
            ENDDO
!---  >          generate p
            DO i = 1,n
               s = 0.0
               DO j = 1,n
                  s = s+f(j)*h(j,i)
               ENDDO
               p(i) = s
            ENDDO
         ENDIF
      ENDIF
      
!-------------------------------------------------------------
!---  >       put xold and y in the correct form for the next step
      DO i = 1,n
         xold(i) = xnew(i)
         y(i) = f(i)
      ENDDO
      
!---  >    if displacements are all less than epsdisp, then converged
      dispmax = 0.0
      DO na = 1,ntype
         nn = 3*(na-1)
         d2 = p(nn+1)**2 + p(nn+2)**2 + p(nn+3)**2
         dispmax = MAX( dispmax, d2)
      ENDDO
      dispmax = xa*SQRT(dispmax)
      IF (dispmax<epsdisp) THEN
         lconv = .TRUE.
      ELSE
         lconv = .FALSE.
      ENDIF
      
!---  >    get new displacements
      DO i = 1,n
         xnew(i) = xold(i)+p(i)
      ENDDO
      DO na = 1,ntype
         nn = 3*(na-1)
         tau0(1,na) = xnew(nn+1)
         tau0(2,na) = xnew(nn+2)
         tau0(3,na) = xnew(nn+3)
      ENDDO
      
      WRITE (6,'(1x/)')
      WRITE (6,*) 'changes in p for step',ist
      DO na = 1,ntype
         nn = 3*(na-1)
         WRITE (6,'(i5,6f12.6)') na,(p(nn+i),i = 1,3),(xnew(nn+j),j=1,3)
      ENDDO
      
      END SUBROUTINE bfgs
      END MODULE m_bfgs
