      MODULE m_boxdim
      CONTAINS
      SUBROUTINE boxdim(
     >                  bmat,
     <                  arltv1,arltv2,arltv3)
c*********************************************************************
c      This subroutine determines the maximum number L, M and N
c      nrgvl, nrgvm, nrgvn, of reciprocal lattice vectors G(L,M,N) 
c      along the directions G(1), G(2), G(3), respectively, for which 
c
c                     | G(L,M,N) | <= GMAX.
c
c      This equation defines a "sphere" and the nrgvl,m,n define
c      the DIMension of the BOX in which the sphere is placed.
c
c      In reality the G(i)'s, do not form a carteasian coordinate 
c      system. Therefore the "sphere" is not a sphere, but an 
c      ellipsoid. For this ellipsoid the largest L, M and N component 
c      is determined as arltv1, arltv2, arltv3 to construct the boxes and
c
c                    nrgvl  = int( gmax/arltv1 ) + 1
c                    nrgvm  = int( gmax/arltv2 ) + 1
c                    nrgvn  = int( gmax/arltv3 ) + 1
c
c      G(i,xyz) is stored in bmat(i,xyz)
c
c      routine by s.bluegel from carpar-program
c 
c                         S. Bl"ugel, IFF, 13. Nov. 97    
c               tested by S. Heinze , IFF, 
c*********************************************************************
      use m_juDFT
c
C     .. Parameters ..
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      REAL,    INTENT (OUT) :: arltv1,arltv2,arltv3
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN)  :: bmat(3,3)
C     ..
C     .. Local Scalars ..
      INTEGER ixyz,j,k
      REAL    denom,eps,one,zero
C     ..
C     .. Local Arrays ..
      REAL det(3,3),rr(3,3)
c     ..
      DATA one,eps,zero/1.0,1e-10,0.0/      
c
c--->  build up quadratic form for ellipsoid
c
      DO j = 1 , 3
         DO k = 1 , 3
            rr(k,j) = zero
            DO ixyz = 1 , 3
               rr(k,j) = rr(k,j) + bmat(k,ixyz)*bmat(j,ixyz)
            ENDDO
         ENDDO
      ENDDO
c
c---> build determinants for Cramer's rule
c
      det(1,1) = rr(2,2)*rr(3,3) - rr(3,2)*rr(2,3)
      det(1,2) = rr(2,1)*rr(3,3) - rr(3,1)*rr(2,3)
      det(1,3) = rr(2,1)*rr(3,2) - rr(3,1)*rr(2,2)
      det(2,1) = rr(1,2)*rr(3,3) - rr(3,2)*rr(1,3)
      det(2,2) = rr(1,1)*rr(3,3) - rr(3,1)*rr(1,3)
      det(2,3) = rr(1,1)*rr(3,2) - rr(3,1)*rr(1,2)
      det(3,1) = rr(1,2)*rr(2,3) - rr(2,2)*rr(1,3)
      det(3,2) = rr(1,1)*rr(2,3) - rr(2,1)*rr(1,3)
      det(3,3) = rr(1,1)*rr(2,2) - rr(2,1)*rr(1,2)
c
c---> check on the zeros of some determinants
c
      DO j = 1 , 3
         IF ( det(j,j) .lt. eps ) THEN
            WRITE (16,
     +                '('' problem with det('',i1,'','',i1,'')'')') j,j
            CALL juDFT_error(" boxdim: determinant",calledby ="boxdim")
         END IF
      ENDDO    
c
c---> scale determinants
c
      DO k = 1 , 3
         denom = one / det(k,k)
         DO j = 1 , 3
            det(k,j) = denom * det(k,j)
         ENDDO
         det(k,k) = one / denom
      ENDDO
c
c---> calculate the maximum l, m, n components of the ellipsoid
c
      arltv1 = sqrt( rr(1,1) + rr(2,2)*det(1,2)*det(1,2)
     >                       + rr(3,3)*det(1,3)*det(1,3)
     >                       - (rr(1,2)+rr(1,2))*det(1,2)
     >                       + (rr(1,3)+rr(1,3))*det(1,3)
     >                       - (rr(2,3)+rr(2,3))*det(1,2)*det(1,3))
      arltv2 = sqrt( rr(2,2) + rr(1,1)*det(2,1)*det(2,1)
     >                       + rr(3,3)*det(2,3)*det(2,3)
     >                       - (rr(1,2)+rr(1,2))*det(2,1)
     >                       + (rr(1,3)+rr(1,3))*det(2,1)*det(2,3)
     >                       - (rr(2,3)+rr(2,3))*det(2,3))
      arltv3 = sqrt( rr(3,3) + rr(1,1)*det(3,1)*det(3,1)
     >                       + rr(2,2)*det(3,2)*det(3,2)
     >                       - (rr(1,2)+rr(1,2))*det(3,1)*det(3,2)
     >                       + (rr(1,3)+rr(1,3))*det(3,1)
     >                       - (rr(2,3)+rr(2,3))*det(3,2))

      RETURN
      END SUBROUTINE
      END
