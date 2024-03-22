module m_sphbessel_integral
contains
   FUNCTION sphbessel_integral(atoms, itype, qnrm, nqnrm, iqnrm1, iqnrm2, l, hybinp, &
                               sphbes0, l_warnin, l_warnout)

      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_atoms), INTENT(IN)   :: atoms

      INTEGER, INTENT(IN)  :: itype, nqnrm, iqnrm1, iqnrm2, l
      REAL, INTENT(IN)  :: qnrm(nqnrm), sphbes0(-1:hybinp%lexp + 2, atoms%ntype, nqnrm)
      LOGICAL, INTENT(IN), OPTIONAL  ::  l_warnin
      LOGICAL, INTENT(INOUT), OPTIONAL  ::  l_warnout
      REAL                  :: sphbessel_integral
      REAL                  :: q1, q2, dq, s, sb01, sb11, sb21, sb31, sb02, sb12
      REAL                  :: sb22, sb32, a1, a2, da, b1, b2, db, c1, c2, dc, r1, r2
      LOGICAL               :: l_warn, l_warned

      IF (PRESENT(l_warnin)) THEN
         l_warn = l_warnin
      ELSE
         l_warn = .TRUE.
      END IF
      l_warned = .FALSE.

      q1 = qnrm(iqnrm1)
      q2 = qnrm(iqnrm2)
      s = atoms%rmt(itype)
      IF (abs(q1) < 1e-12 .AND. abs(q2) < 1e-12) THEN
      IF (l > 0) THEN
         sphbessel_integral = 0
      ELSE
         sphbessel_integral = 2*s**5/15
      ENDIF
      ELSE IF (abs(q1) < 1e-12 .OR. abs(q2) < 1e-12) THEN
      IF (l > 0) THEN
         sphbessel_integral = 0
      ELSE IF (abs(q1) < 1e-12) THEN
         sphbessel_integral = s**3/(3*q2**2)*(q2*s*sphbes0(1, itype, iqnrm2) &
                                              + sphbes0(2, itype, iqnrm2))
      ELSE
         sphbessel_integral = s**3/(3*q1**2)*(q1*s*sphbes0(1, itype, iqnrm1) &
                                              + sphbes0(2, itype, iqnrm1))
      ENDIF
      ELSE IF (abs(q1 - q2) < 1e-12) THEN
      sphbessel_integral = s**3/(2*q1**2)*((2*l + 3)*sphbes0(l + 1, itype, iqnrm1)**2 - &
                                           (2*l + 1)*sphbes0(l, itype, iqnrm1)*sphbes0(l + 2, itype, iqnrm1))
      ELSE ! We use either if two fromulas that are stable for high and small q1/q2 respectively
      sb01 = sphbes0(l - 1, itype, iqnrm1)
      sb11 = sphbes0(l, itype, iqnrm1)
      sb21 = sphbes0(l + 1, itype, iqnrm1)
      sb31 = sphbes0(l + 2, itype, iqnrm1)
      sb02 = sphbes0(l - 1, itype, iqnrm2)
      sb12 = sphbes0(l, itype, iqnrm2)
      sb22 = sphbes0(l + 1, itype, iqnrm2)
      sb32 = sphbes0(l + 2, itype, iqnrm2)
      dq = q1**2 - q2**2
      a1 = q2/q1*sb21*sb02
      a2 = q1/q2*sb22*sb01
      da = a1 - a2
      b1 = sb31*sb12
      b2 = sb32*sb11
      db = b1 - b2
      c1 = sb21*sb22/(q1*q2)
      c2 = db/dq*(2*l + 1)/(2*l + 3)
      dc = c1 + c2
      r1 = ABS(da/a1)
      r2 = MIN(ABS(db/b1), ABS(dc/c1))
! Ensure numerical stability. If both formulas are not sufficiently stable, the program stops.
      IF (r1 > r2) THEN
      IF (r1 < 1e-6 .AND. l_warn) THEN
         WRITE (oUnit, '(A,E12.5,A,E12.5,A)') 'sphbessel_integral: Warning! Formula One possibly unstable. Ratios:', &
            r1, '(', r2, ')'
         WRITE (oUnit, '(A,2F15.10,I4)') '                    Current qnorms and atom type:', q1, q2, itype
         l_warned = .TRUE.
      END IF
      sphbessel_integral = s**3/dq*da
      ELSE
      IF (r2 < 1e-6 .AND. l_warn) THEN
         WRITE (oUnit, '(A,E13.5,A,E13.5,A)') 'sphbessel_integral: Warning! Formula Two possibly unstable. Ratios:', &
            r2, '(', r1, ')'
         WRITE (oUnit, '(A,2F15.10,I4)') '                    Current qnorms and atom type:', &
            q1, q2, itype
         l_warned = .TRUE.
      END IF
      sphbessel_integral = s**3*dc
      END IF
      END IF

      IF (PRESENT(l_warnout)) l_warnout = l_warned

   END FUNCTION sphbessel_integral
end module m_sphbessel_integral
