      MODULE m_grdrsvac
      use m_juDFT
c----------------------------------------------------------------------
c     for a function defined on a vacuum layer
c     the in-plane derivatives are evaluated in real space
c
c     based on 'rhzgrd' coded by t.asada. june,1995.
c---------------------------------------------------------------------
      CONTAINS
      SUBROUTINE grdrsvac( 
     >                    ro,bmat,xmax1,xmax2,ndvgrd,
     <                    dxro,dyro)

c----------------------------------------------------------------------
c
c input: 
c  ro(0:xmax1*xmax2-1)  
c   any quantity stored in the 2-dim box (xmax1 x xmax2)
c  bmat
c   bravais matrix of reciprocal space 
c  ndvgrd
c   number of points used when calculating derivative (3 <= ndvgrd <= 6)
c       
c output:
c  dxro(0:xmax1*xmax2-1) , dxro(0:xmax1*xmax2-1) 
c   (d ro / d x) , (d ro /d y)  in non-internal coordinates 
c       
c----------------------------------------------------------------------
c
      USE m_constants, ONLY : pimach
      IMPLICIT NONE 
c     ..
c     .. Scalar arguments ..
      INTEGER, INTENT (IN) :: xmax1,xmax2,ndvgrd
c     ..
c     .. Array arguments ..
      REAL,    INTENT (IN)  :: ro(0:xmax1*xmax2-1)
      REAL,    INTENT (IN)  :: bmat(3,3) 
c     ..
c     .. Array output ..
      REAL,    INTENT (OUT) :: dxro(0:xmax1*xmax2-1)
      REAL,    INTENT (OUT) :: dyro(0:xmax1*xmax2-1)   
c     ..
c     .. Locals ..
      INTEGER :: xmax(2)
      INTEGER :: direction,xyz(2),x1,x2,i,ii(-3:2) 
      REAL    :: pi  
      REAL    :: dx
      REAL    :: drointern(0:xmax1*xmax2-1,2)  

      pi = pimach()

      xmax(1)= xmax1
      xmax(2)= xmax2
      DO i=1,2
        IF ( xmax(i) < 3 ) THEN
           CALL juDFT_error("grid to small",calledby="grdrsvac")
        END IF
      END DO 
      IF ( (ndvgrd < 3) .or. (ndvgrd > 6) ) THEN
         CALL juDFT_error("ndvgrd notin [3,6]",calledby="grdrsvac")
      ENDIF


      DO direction=1,2 

        dx= 1./REAL(xmax(direction)) 

        DO x1=0,xmax(1)-1
          DO x2=0,xmax(2)-1

            DO i= -3,2
              xyz(1)= x1
              xyz(2)= x2
              xyz(direction)= xyz(direction)+i 
              ! make use of periodic boundary cond. in interstitial: 
              IF ( xyz(direction) < 0 ) THEN
                xyz(direction)= xyz(direction)+xmax(direction)
              END IF 
              IF ( xyz(direction) >= xmax(direction) ) THEN
                xyz(direction)= xyz(direction)-xmax(direction) 
              END IF 
              ! find coordinates in 1-dim array ro:
              ii(i)= xyz(2)*xmax(1) + xyz(1) 
            END DO

            IF (ndvgrd.EQ.3) THEN
              drointern(ii(0),direction)=  
     &          df3( ro(ii(-1)), 
     &               ro(ii(0)),ro(ii(1)), dx)
            ELSEIF (ndvgrd.EQ.4) THEN
              drointern(ii(0),direction)= 
     &          df4( ro(ii(-1)),
     &               ro(ii(0)),ro(ii(1)),ro(ii(2)), dx)
            ELSEIF (ndvgrd.EQ.5) THEN
              drointern(ii(0),direction)= 
     &           df5( ro(ii(-2)),ro(ii(-1)),
     &                ro(ii(0)),ro(ii(1)),ro(ii(2)), dx)
            ELSEIF (ndvgrd.EQ.6) THEN
              drointern(ii(0),direction)= 
     &          df6( ro(ii(-3)),ro(ii(-2)),ro(ii(-1)),
     &               ro(ii(0)),ro(ii(1)),ro(ii(2)), dx)
            ENDIF

          END DO
        END DO 

      END DO

 
      DO i=0,xmax(1)*xmax(2)-1 

        dxro(i)=   bmat(1,1)*drointern(i,1) 
     &           + bmat(2,1)*drointern(i,2)
        dxro(i)= dxro(i)/(2.*pi) 
        dyro(i)=   bmat(1,2)*drointern(i,1) 
     &           + bmat(2,2)*drointern(i,2)
        dyro(i)= dyro(i)/(2.*pi) 

      END DO     

      END SUBROUTINE grdrsvac
!--------------------------------------------------------------------
! Functions: formulae for 1st deriv.:
!
      REAL FUNCTION df3(g1,f0,f1,d)             ! three point formula
        REAL g1,f0,f1,d
        df3 = (-1*g1-0*f0+f1)/ (2*d)
      END FUNCTION df3

      REAL FUNCTION df4(g1,f0,f1,f2,d)          ! four point formula
        REAL g1,f0,f1,f2,d
        df4 = (-2*g1-3*f0+6*f1-f2)/ (6*d)
      END FUNCTION df4

      REAL FUNCTION df5(g2,g1,f0,f1,f2,d)       ! five point formula
        REAL g2,g1,f0,f1,f2,d
        df5 = (2*g2-16*g1-0*f0+16*f1-2*f2)/ (24*d)
      END FUNCTION df5

      REAL FUNCTION df6(g3,g2,g1,f0,f1,f2,d)   ! six point formula
        REAL g3,g2,g1,f0,f1,f2,d
        df6 = (-4*g3+30*g2-120*g1+40*f0+60*f1-6*f2)/ (120*d)
      END FUNCTION df6

!----------------------------------------------------------------------
      END MODULE m_grdrsvac

