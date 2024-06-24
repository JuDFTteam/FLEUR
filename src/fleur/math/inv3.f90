MODULE m_inv3
  !-----------------------------------
  !     invert 3x3 matrix
  !-----------------------------------

  PRIVATE

  INTERFACE inv3
     MODULE PROCEDURE inv3i,inv3r
  END INTERFACE inv3

  PUBLIC :: inv3

CONTAINS
  SUBROUTINE inv3r(a,b,d)

    IMPLICIT NONE
    !     ..
    !     .. Arguments ..
    REAL, INTENT (IN)  :: a(3,3)
    REAL, INTENT (OUT) :: b(3,3)  ! inverse matrix
    REAL, INTENT (OUT) :: d       ! determinant
    !     ..
    d = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + &
         a(2,1)*a(3,2)*a(1,3) - a(1,3)*a(2,2)*a(3,1) - &
         a(2,3)*a(3,2)*a(1,1) - a(2,1)*a(1,2)*a(3,3)
    b(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))/d
    b(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))/d
    b(1,3) = (a(1,2)*a(2,3)-a(2,2)*a(1,3))/d
    b(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))/d
    b(2,2) = (a(1,1)*a(3,3)-a(3,1)*a(1,3))/d
    b(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))/d
    b(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))/d
    b(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))/d
    b(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/d

  END SUBROUTINE inv3r

  SUBROUTINE inv3i(a,b,d)

    IMPLICIT NONE
    !     ..
    !     .. Arguments ..
    INTEGER, INTENT (IN)  :: a(3,3)
    INTEGER, INTENT (OUT) :: b(3,3)  ! inverse matrix
    INTEGER, INTENT (OUT) :: d       ! determinant
    !     ..
    d = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) +&
         a(2,1)*a(3,2)*a(1,3) - a(1,3)*a(2,2)*a(3,1) -&
         a(2,3)*a(3,2)*a(1,1) - a(2,1)*a(1,2)*a(3,3)
    b(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))/d
    b(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))/d
    b(1,3) = (a(1,2)*a(2,3)-a(2,2)*a(1,3))/d
    b(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))/d
    b(2,2) = (a(1,1)*a(3,3)-a(3,1)*a(1,3))/d
    b(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))/d
    b(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))/d
    b(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))/d
    b(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/d

  END SUBROUTINE inv3i
END MODULE m_inv3
