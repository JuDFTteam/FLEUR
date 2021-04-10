MODULE m_angles

   USE m_types_sym
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   !Calculate the correct phase shifts for the rotation in spin-space
   !(needed for LDA+U and Green's functions with noco%l_mperp)

   SUBROUTINE angles(sym)

      TYPE(t_sym),   INTENT(INOUT)  :: sym

      INTEGER iop,d,t

      sym%phase = 0.0

      DO iop = 1, sym%nop
         d = det(sym%mrot(:,:,iop))
         t = (sym%mrot(1,1,iop)+sym%mrot(2,2,iop)+sym%mrot(3,3,iop)) * d
         IF(t.EQ.-1) sym%phase(iop) = 1.0
         IF(t.EQ.0)  sym%phase(iop) = 2.0/3.0
         IF(t.EQ.1)  sym%phase(iop) = 1.0/2.0
         IF(t.EQ.2)  sym%phase(iop) = 1.0/3.0
         IF(t.EQ.3)  sym%phase(iop) = 0.0
         sym%phase(iop) = d*sym%phase(iop)*pi_const
      ENDDO

   END SUBROUTINE

   INTEGER FUNCTION det(m)
      INTEGER m(:,:)
      det = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
            m(2,1)*m(3,2)*m(1,3) - m(1,3)*m(2,2)*m(3,1) - &
            m(2,3)*m(3,2)*m(1,1) - m(2,1)*m(1,2)*m(3,3)
   END FUNCTION det


END MODULE m_angles