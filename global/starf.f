      MODULE m_starf
!     ********************************************************
!     synthesize 2-d and 3-d star function at point r, which
!     is given in internal coordinates
!     ********************************************************    
      USE m_constants
      USE m_spgrot
      IMPLICIT NONE
      CONTAINS
!     ********************************************************
      SUBROUTINE starf2(
     >                  nop2,ng2,kv2,mrot,symor,tau,r,invtab,
     <                  sf)
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nop2,ng2
      LOGICAL, INTENT (IN) :: symor
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: kv2(2,ng2),mrot(3,3,nop2)
      INTEGER, INTENT (IN) :: invtab(nop2)
      REAL,    INTENT (IN) :: r(3),tau(3,nop2)
      COMPLEX, INTENT (OUT):: sf(ng2)
!     ..
!     .. Local Arrays ..
      INTEGER kr(3,nop2),kv(3),k,n
      REAL    arg
      COMPLEX ph(nop2)

    
      DO k = 1,ng2
         kv(1) = kv2(1,k)
         kv(2) = kv2(2,k)
         kv(3) = 0
         sf(k) = 0.0

         CALL spgrot(
     >               nop2,symor,mrot,tau,invtab,
     >               kv,
     <               kr,ph)

         DO n = 1,nop2
            arg = tpi_const* (kr(1,n)*r(1)+kr(2,n)*r(2))
            sf(k) = sf(k) + ph(n) * cmplx(cos(arg),sin(arg))
         ENDDO
         sf(k) = sf(k)/nop2
      ENDDO

      END SUBROUTINE starf2
!     ********************************************************
      SUBROUTINE starf3(
     >                  nop,ng3,symor,kv3,mrot,tau,r,invtab,
     <                  sf)
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nop,ng3
      LOGICAL, INTENT (IN) :: symor
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: kv3(3,ng3),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invtab(nop)

      REAL,    INTENT (IN) :: tau(3,nop),r(3)
      COMPLEX, INTENT (OUT):: sf(ng3)
!     ..
!     .. Local Arrays ..
      INTEGER kr(3,nop),k,n
      REAL    arg
      COMPLEX ph(nop)

    
      DO k = 1,ng3
         CALL spgrot(
     >               nop,symor,mrot,tau,invtab,
     >               kv3(1,k),
     <               kr,ph)
         sf(k) = 0.0
         DO n = 1,nop
            arg = tpi_const* (kr(1,n)*r(1)+kr(2,n)*r(2)+kr(3,n)*r(3))
            sf(k) = sf(k) + ph(n) * cmplx(cos(arg),sin(arg))
         ENDDO
         sf(k) = sf(k)/nop
      ENDDO

      END SUBROUTINE starf3
!     ********************************************************
      END MODULE m_starf
