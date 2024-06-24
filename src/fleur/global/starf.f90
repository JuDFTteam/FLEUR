MODULE m_starf
    ! Construct the 2D and 3D star functions for a real space point r
    ! given in internal coordinates.
    ! Formula:
    ! sf(k) = 1/N_{op} \sum_{op(k)} e^{i * (\mat{R}_{op}\vec{G}_{k})^T
    !                                      \cdot (\vec{r} - \vec{\tau}_{op})}

    USE m_constants
    USE m_spgrot

    IMPLICIT NONE

CONTAINS
    SUBROUTINE starf2(nop2, ng2, kv2, mrot, symor, tau, r, invtab, sf, center)

        ! Scalar arguments:
        INTEGER, INTENT(IN) :: nop2, ng2
        LOGICAL, INTENT(IN) :: symor

        ! Array arguments:
        INTEGER, INTENT(IN)  :: kv2(2, ng2), mrot(3, 3, nop2)
        INTEGER, INTENT(IN)  :: invtab(nop2)

        REAL,    INTENT(IN)  :: r(3), tau(3, nop2)
        COMPLEX, INTENT(OUT) :: sf(ng2)
        REAL, INTENT(IN), OPTIONAL :: center(3)

        ! Local scalars:
        INTEGER :: k, n
        REAL    :: arg

        ! Local arrays:
        INTEGER :: kr(3, nop2), kv(3)
        COMPLEX :: ph(nop2)

        DO k = 1,ng2
            kv(1) = kv2(1, k)
            kv(2) = kv2(2, k)
            kv(3) = 0

            CALL spgrot(nop2, symor, mrot, tau, invtab, kv, kr, ph)

            sf(k) = 0.0

            DO n = 1, nop2
                IF (.NOT. PRESENT(center)) THEN
                    arg = tpi_const * ( (kr(1, n))* r(1) + kr(2, n) * r(2) ) 
                ELSE
                    arg = tpi_const * ( (kr(1, n) + center(1) )* r(1) + (kr(2, n) + center(2)) * r(2) ) 
                END IF 
                ! Sum up e^{i * \vec{G}_{op}^T
                !               \cdot (\vec{r} - \vec{\tau}_{op})}
                    sf(k) = sf(k) + ph(n) * cmplx( cos(arg), sin(arg) )
            END DO

            sf(k) = sf(k) / nop2

        END DO

    END SUBROUTINE starf2

    SUBROUTINE starf3(nop,ng3,symor,kv3,mrot,tau,r,invtab,sf,center)

        ! Scalar arguments:
        INTEGER, INTENT(IN) :: nop, ng3
        LOGICAL, INTENT(IN) :: symor

        ! Array arguments:
        INTEGER, INTENT(IN) :: kv3(3, ng3), mrot(3, 3, nop)
        INTEGER, INTENT(IN) :: invtab(nop)

        REAL,    INTENT(IN)   :: tau(3, nop), r(3)
        COMPLEX, INTENT (OUT) :: sf(ng3)
        REAL, INTENT(IN), OPTIONAL :: center(3)
        ! Local scalars:
        INTEGER :: k,n
        REAL    :: arg

        ! Local arrays:
        INTEGER :: kr(3,nop)
        COMPLEX :: ph(nop)

        DO k = 1, ng3

            CALL spgrot(nop, symor, mrot, tau, invtab, kv3(:, k), kr, ph)

            sf(k) = 0.0

            DO n = 1, nop
                IF (PRESENT(center)) THEN
                    arg = tpi_const * dot_product( real(kr(:, n)) + center , r )   !! if flags + center
                ELSE
                    arg = tpi_const * dot_product( real(kr(:, n)) , r ) 
                END IF 
                ! Sum up e^{i * \vec{G}_{op}^T
                !               \cdot (\vec{r} - \vec{\tau}_{op})}
                sf(k) = sf(k) + ph(n) * cmplx( cos(arg), sin(arg) )   
            END DO

            sf(k) = sf(k) / nop

        END DO

    END SUBROUTINE starf3

END MODULE m_starf
