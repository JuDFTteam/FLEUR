      MODULE m_od_chirot
      CONTAINS
      SUBROUTINE od_chirot(
     >                  odi,ods,bmat,k,
     <                  kr,phas)

****************************************************************
*     perform symmetry operations of the cylindrical group 
*     of symmetries on the reciprocal vector k, given in internal
*     coordinates
*                                       Y.Mokrousov
****************************************************************

      USE m_types, ONLY : od_inp, od_sym

      IMPLICIT NONE
!     ..
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
      INTEGER, INTENT (IN)  :: k(3)
      REAL,    INTENT (IN)  :: bmat(3,3)
      REAL,    INTENT (OUT) :: kr(3,ods%nop)
      COMPLEX, INTENT (OUT) :: phas(ods%nop)
!     ..
!     .. Local Scalars ..
      INTEGER n,inv
!     ..
!     ..

      DO n = 1,ods%nop
         kr(1,n) = k(1)*ods%mrot(1,1,n) +
     +             k(2)*ods%mrot(2,1,n) +
     +             k(3)*ods%mrot(3,1,n)
         kr(2,n) = k(1)*ods%mrot(1,2,n) + 
     +             k(2)*ods%mrot(2,2,n) +
     +             k(3)*ods%mrot(3,2,n)
         kr(3,n) = k(1)*ods%mrot(1,3,n) + 
     +             k(2)*ods%mrot(2,3,n) +
     +             k(3)*ods%mrot(3,3,n)
      ENDDO

      IF (odi%chi.NE.1) THEN
         DO n = 1,ods%nop
            inv = ods%invtab(n)
            phas(n) = exp(cmplx(0.,bmat(3,3)*kr(3,n)*ods%tau(3,inv)))
         END DO
      ELSE
         DO n = 1,ods%nop
            phas(n) = cmplx(1.,0.)
         END DO
      END IF
      
      RETURN
      END SUBROUTINE od_chirot
      END MODULE m_od_chirot

 
