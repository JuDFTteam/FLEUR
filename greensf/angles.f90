MODULE m_angles

CONTAINS

   !Calculate the correct angles for the rotation in spin-space (needed for LDA+U and Greens functions with noco)

   SUBROUTINE angles(sym,angle)

      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_sym),   INTENT(IN)  :: sym
      REAL,          INTENT(OUT) :: angle(:)

      INTEGER iop,d,t

      angle = 0.0

      DO iop = 1, sym%nop
         d = det(sym%mrot(:,:,iop))
         t = (sym%mrot(1,1,iop)+sym%mrot(2,2,iop)+sym%mrot(3,3,iop)) * d
         IF(t.EQ.-1) angle(iop) = 1.0
         IF(t.EQ.0)  angle(iop) = 2.0/3.0
         IF(t.EQ.1)  angle(iop) = 1.0/2.0
         IF(t.EQ.2)  angle(iop) = 1.0/3.0
         IF(t.EQ.3)  angle(iop) = 0.0
         angle(iop) = d*angle(iop)*pi_const
      ENDDO

   END SUBROUTINE

   INTEGER FUNCTION det(m)
      INTEGER m(:,:)
      det = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
            m(2,1)*m(3,2)*m(1,3) - m(1,3)*m(2,2)*m(3,1) - &
            m(2,3)*m(3,2)*m(1,1) - m(2,1)*m(1,2)*m(3,3)
   END FUNCTION det


END MODULE m_angles