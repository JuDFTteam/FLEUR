module m_intgrf
   implicit none
   TYPE :: intgrf_out
      REAL      :: value    ! value of the integration
      INTEGER :: ierror   ! error code
   END TYPE intgrf_out

   !     error and warning codes for intgrf function
   INTEGER, PARAMETER :: NO_ERROR = 0
   INTEGER, PARAMETER :: NEGATIVE_EXPONENT_WARNING = 1
   INTEGER, PARAMETER :: NEGATIVE_EXPONENT_ERROR = 2
CONTAINS

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   ! Integrates function f numerically (Lagrange and Simpson integration)
   ! on grid(itype) and is much faster than intgr.
   ! (Only normal outward integration.)
   ! Before first use of this function it has to be initialized with
   ! intgrf_init.

   FUNCTION intgrf(f, atoms, itype, gridf)
      use m_juDFT
      use m_types_setup
      IMPLICIT NONE

      REAL                 :: intgrf
      type(t_atoms)        :: atoms
      INTEGER, INTENT(IN)  :: itype
      REAL, INTENT(IN)  :: gridf(atoms%jmtd, atoms%ntype)
      REAL, INTENT(IN)  :: f(:)
      !     - local -
      TYPE(intgrf_out)     :: fct_res

      fct_res = pure_intgrf(f, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, &
                            atoms%ntype, itype, gridf)
      IF (fct_res%ierror == NEGATIVE_EXPONENT_WARNING) THEN
         write (6, *) 'intgrf: Warning!'// &
            'Negative exponent x in extrapolation a+c*r**x'
      ELSEIF (fct_res%ierror == NEGATIVE_EXPONENT_ERROR) THEN
         write (6, *) &
            'intgrf: Negative exponent x in extrapolation a+c*r**x'
         CALL juDFT_error( &
            'intgrf: Negative exponent x in extrapolation a+c*r**x')
      END IF
      intgrf = fct_res%value

   END FUNCTION intgrf

   !     pure wrapper for intgrf with same functionality
   !     can be used within forall loops

   PURE FUNCTION pure_intgrf(f, jri, jmtd, rmsh, dx, ntype, itype, gridf)

      IMPLICIT NONE

      TYPE(intgrf_out)     :: pure_intgrf

      INTEGER, INTENT(IN)  :: itype, ntype, jmtd
      INTEGER, INTENT(IN)  :: jri(ntype)
      REAL, INTENT(IN)  :: dx(ntype), rmsh(jmtd, ntype)
      REAL, INTENT(IN)  :: gridf(jmtd, ntype)
      REAL, INTENT(IN)  :: f(:)
      !     - local -
      INTEGER                :: n
      REAL                   :: r1, h, a, x

      n = jri(itype)
      r1 = rmsh(1, itype)
      h = dx(itype)

      pure_intgrf%ierror = NO_ERROR

      ! integral from 0 to r1 approximated by leading term in power series expansion
      IF (f(1)*f(2) > 1e-10 .AND. h > 0) THEN
         IF (f(2) == f(1)) THEN
            pure_intgrf%value = r1*f(1)
         ELSE
            x = (f(3) - f(2))/(f(2) - f(1))
            a = (f(2) - x*f(1))/(1 - x)
            x = log(x)/h
            IF (x < 0) THEN
               IF (x > -1) THEN
                  pure_intgrf%ierror = NEGATIVE_EXPONENT_WARNING
               ELSE IF (x <= -1) THEN
                  pure_intgrf%ierror = NEGATIVE_EXPONENT_ERROR
                  RETURN
               END IF
            END IF

            pure_intgrf%value = r1*(f(1) + x*a)/(x + 1)

            !           x      = f(2) / f(1)
            !           x      = log(x)/h
            !           IF(x.lt.0) THEN
            !             IF(x.gt.-1) write(6,'(A,ES9.1)') 'intgrf: Warning!&
            !      &                                        Negative exponent x in&
            !      &                                        extrapolation c*r**x:',x
            !             IF(x.le.-1) write(6,'(A,ES9.1)') 'intgrf: Negative exponent&
            !      &                                        x in extrapolation&
            !      &                                        c*r**x:',x
            !             IF(x.le.-1) call juDFT_error('intgrf: Negative exponent&
            !      &                                        x in extrapolation &
            !      &                                        c*r**x')
            !           END IF
            !           intgrf = (r1*f(1))/(x+1)

         END IF
      ELSE
         pure_intgrf%value = 0
      END IF

      ! integrate from r(1) to r(n) by multiplying with gridf
      pure_intgrf%value = pure_intgrf%value &
                          + dot_product(gridf(:n, itype), f(:n))

   END FUNCTION pure_intgrf

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   !     Initializes fast numerical integration intgrf (see below).

   SUBROUTINE intgrf_init(ntype, jmtd, jri, dx, rmsh, gridf)

      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: ntype, jmtd
      INTEGER, INTENT(IN)   :: jri(ntype)
      REAL, INTENT(IN)   :: dx(ntype), rmsh(jmtd, ntype)
      REAL, INTENT(OUT), ALLOCATABLE  :: gridf(:, :)

      !     - local -
      INTEGER                :: i, j, itype
      INTEGER                :: n, nstep, n0 = 6
      INTEGER, PARAMETER   :: simpson(7) = (/41, 216, 27, 272, 27, 216, 41/)
      REAL                   :: r1, h, dr
      REAL                   :: r(7)
      REAL, PARAMETER      :: lagrange(7, 6) = reshape( &
      (/19087., 65112., -46461., 37504., -20211., 6312., -863., &
      -863., 25128., 46989., -16256., 7299., -2088., 271., &
      &               271., -2760., 30819., 37504., -6771., 1608., -191., &
      -191., 1608., -6771., 37504., 30819., -2760., 271., &
      &               271., -2088., 7299., -16256., 46989., 25128., -863., &
      -863., 6312., -20211., 37504., -46461., 65112., 19087./), & ! The last row is actually never used.
      (/7, 6/))

      n = jmtd
      ALLOCATE (gridf(n, ntype))

      gridf = 0

      DO itype = 1, ntype

         n = jri(itype)
         r1 = rmsh(1, itype)
         h = dx(itype)

         nstep = (n - 1)/6
         n0 = n - 6*nstep
         dr = exp(h)

         ! Calculate Lagrange-integration coefficients from r(1) to r(n0)
         r(1) = r1
         IF (n0 > 1) THEN
            DO i = 2, 7
               r(i) = r(i - 1)*dr
            END DO
            DO i = 1, 7
               gridf(i, itype) = h/60480*r(i)*sum(lagrange(i, 1:n0 - 1))
            END DO
            r(1) = r(n0)
         END IF

         ! Calculate Simpson-integration coefficients from r(n0) to r(n)
         DO i = 1, nstep
            DO j = 2, 7
               r(j) = r(j - 1)*dr
            END DO
            DO j = n0, n0 + 6
               gridf(j, itype) = gridf(j, itype) + h/140*r(j - n0 + 1)* &
                                 simpson(j - n0 + 1)
            END DO
            n0 = n0 + 6
            r(1) = r(7)
         END DO

      END DO

   END SUBROUTINE intgrf_init
end module m_intgrf
