      MODULE m_spgrot
!     **********************************************************
!     perform space group operations of film
!     **********************************************************
      CONTAINS
      SUBROUTINE spgrot(
     >                  nop,symor,mrot,tau,invtab,
     >                  k,
     <                  kr,phas)
      USE m_constants
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)  :: nop
      LOGICAL, INTENT (IN)  :: symor
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: k(3),mrot(3,3,nop),invtab(nop)
      REAL,    INTENT (IN)  :: tau(3,nop)
      INTEGER, INTENT (OUT) :: kr(3,nop)
      COMPLEX,OPTIONAL, INTENT (OUT) :: phas(nop)
!     ..
!     .. Local Scalars ..
      INTEGER n,ni

      DO n = 1,nop
         kr(:,n) = matmul(k,mrot(:,:,n))
      ENDDO
      if (.not. present(phas)) RETURN
      IF (symor) THEN
            phas(:) = 1.
      ELSE
         DO n = 1,nop
            ni = invtab(n)
            phas(n) = exp(cmplx(0.0,-1.0)*tpi_const* 
     +                   dot_product(real(kr(:,n)),tau(:,ni)))
! note that, in general phas(n) could be complex!
         ENDDO
      END IF
      RETURN
      END SUBROUTINE spgrot
      END MODULE m_spgrot
