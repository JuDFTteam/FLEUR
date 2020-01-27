MODULE m_bfgs0

CONTAINS

!*******************************************************************
!     checks whether a file exists with positions and forces to
!     use to relax structure using the BFGS quasi-Newton method.
!     istep0 is the number of steps included on file.
!     h,y, and xold will be saved correctly for further use.
!*******************************************************************
SUBROUTINE bfgs0(ntype,istep0,xold,y,h)
      
   IMPLICIT NONE

   ! Scalar Arguments
   INTEGER, INTENT (IN) :: ntype
   INTEGER, INTENT (OUT):: istep0

   ! Array Arguments
   REAL,    INTENT (OUT):: y(3*ntype),xold(3*ntype)
   REAL,    INTENT (OUT):: h(3*ntype,3*ntype)

   ! Local Scalars
   INTEGER n,is,i,j
   REAL yy,py,s,gamma
   CHARACTER :: testchar
   INTEGER   :: maxhist
   LOGICAL   :: l_skip

   ! Local Arrays
   REAL f(3*ntype),p(3*ntype),v(3*ntype),xnew(3*ntype)

   ! open the file with forces

   OPEN (43,file='forces.dat',status='unknown',form='formatted')
   REWIND 43

   n=3*ntype
   istep0=0

   ! There might be a line like
   ! maxhist=3
   ! at the beginning of the file
   READ(43,'(a)',end=100) testchar
   REWIND 43
   IF (testchar =='m') THEN
      READ(43,"(8x,i5)") maxhist
      DO
         READ(43,*,END = 99) !read the line containing the energy
         istep0 = istep0+1
         DO j = 1,ntype
            READ(43,*,err= 100, end=100)
         END DO
      END DO
      99 REWIND 43
      READ(43,*) !skip maxhist line
      !skip all but last maxhist entries
      DO i = 1, MAX(istep0-maxhist,0)
         READ(43,*,end=100)
         DO j = 1,ntype
            READ(43,*,err= 100, end=100)
         END DO
      END DO
   END IF
   istep0 = 0

   DO ! read from file until end of data
      l_skip = .FALSE.
      READ(43,'(33x,l1)',END = 101,err = 101) l_skip
      101 READ (43,'(3f16.9,3f14.9)',END = 100,err = 100) ((xnew(i+3*j),i = 1,3),(f(i+3*j),i = 1,3),j = 0,ntype-1) 
      IF (l_skip) CYCLE !skip this entry

      istep0 = istep0+1

      IF (istep0.EQ.1) THEN
         ! first step

         ! set h to identity matrix
         h = 0.0
         DO j = 1,n
            h(j,j) = 1.0
         END DO

      ELSE
         ! k-th step
            
         ! determine p and shift x as needed
         p(:) = xnew(:) - xold(:)
         ! update the change in gradients
         y(:) = y(:) - f(:)
         ! get necessary inner products and H|y>
         py = dot_product(p,y)
         yy = 0.0
         DO i = 1,n
            s = 0.0
            DO j = 1,n
               s = s+y(j)*h(j,i)
            END DO
            v(i) = s
            yy = yy+y(i)*s
         END DO
         ! check that update will leave h positive definite, if not, restart with diagonal matrix
         IF (py<=0.0) THEN
            WRITE (6,*) 'bfgs0: <p|y> < 0, restart istep0 =',istep0
            h = 0.0
            DO j = 1,n
               h(j,j) = 1.0
            END DO
            istep0 = 1
         ELSE
            ! for second step, use oren-spedicato scaling
            ! for standard scaling, use gamma = 1.
            IF (istep0  == 2) THEN
               gamma = py/yy
            ELSE
               gamma = 1.0
            END IF
            ! update h
            DO j = 1,n
               DO i = 1,n
                  h(i,j) = (h(i,j) - (v(i)*p(j)+p(i)*v(j))/py)*gamma + (1.+gamma*yy/py)*p(i)*p(j)/py
               END DO
            END DO
         END IF
      END IF

      ! put y in the correct form for the next step
      DO i = 1,n
         xold(i) = xnew(i)
         y(i) = f(i)
      END DO
         
   END DO
         
   100 CLOSE(43)

   IF (istep0==0) THEN
      WRITE (6,1000)
   ELSE
      WRITE (6,1001) istep0
   END IF

   1000 FORMAT(' bfgs0: No previous atomic relaxation data found'/)
   1001 FORMAT(' bfgs0: Approximate hessian constructed using',i4,' steps'/)

   END SUBROUTINE bfgs0

END MODULE m_bfgs0
