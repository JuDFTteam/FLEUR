MODULE m_spgrot
    ! Perform the space group operations of the system.
    ! I.e. construct all G vectors (and optionally phases) of a star from
    ! its representative reciprocal lattice vector k in internal coordinates.

CONTAINS
    SUBROUTINE spgrot(nop, symor, mrot, tau, invtab, k, kr, phas)

        USE m_constants

        IMPLICIT NONE

        ! Scalar arguments:
        INTEGER, INTENT(IN)  :: nop
        LOGICAL, INTENT(IN)  :: symor

        ! Array arguments:
        INTEGER, INTENT(IN)  :: k(3), mrot(3, 3, nop), invtab(nop)
        REAL,    INTENT(IN)  :: tau(3, nop)

        INTEGER,           INTENT(OUT) :: kr(3, nop)
        COMPLEX, OPTIONAL, INTENT(OUT) :: phas(nop) ! Could be complex!

        ! Local scalars:
        INTEGER :: n, ni

        DO n = 1, nop
            ! Construct \vec{G}_{op}^T = \vec{G}_{star}^T \cdot \mat{R}_{op}
            kr(:, n) = matmul(k, mrot(:, :, n))
        END DO

        IF (.NOT.PRESENT(phas)) RETURN

        IF (symor) THEN
            phas(:) = 1.0
        ELSE
            DO n = 1,nop
                ni = invtab(n)
                ! Construct e^{-i * \vec{G}_{op}^T \cdot \vec{\tau}_{op}}
                phas(n) = exp( -ImagUnit * tpi_const &
                      & * dot_product(real(kr(:, n)), tau(:, ni)) )
            END DO
        END IF

    END SUBROUTINE spgrot
END MODULE m_spgrot
