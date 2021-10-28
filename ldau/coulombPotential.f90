module m_coulombPotential

    use m_types
    use m_uj2f
    use m_umtx

    implicit none

    contains

    subroutine coulombPotential(density, ldau, jspins, l_spinoffd, potential, interaction_energy)

        !This subroutine calculates the DFT+U potential matrix excluding the double counting
    
        COMPLEX,          INTENT(IN)     :: density(-lmaxU_const:,-lmaxU_const:,:)
        TYPE(t_utype),    INTENT(IN)     :: ldau
        INTEGER,          INTENT(IN)     :: jspins
        LOGICAL,          INTENT(IN)     :: l_spinoffd
        COMPLEX,          INTENT(INOUT)  :: potential(-lmaxU_const:,-lmaxU_const:,:)
        REAL,             INTENT(OUT)    :: interaction_energy

        INTEGER :: m,mp,q,p,ispin,jspin,spin_dim
        REAL    :: spin_deg,f0,f2,f4,f6,energy_contribution,total_charge
        COMPLEX :: vu
        REAL,    ALLOCATABLE :: umatrix(:,:,:,:)
        COMPLEX, ALLOCATABLE :: Vdc(:,:,:)

        spin_dim = SIZE(density,3)
        spin_deg = 1.0 / (3 - jspins)

        IF(SIZE(potential,3) /= spin_dim) CALL juDFT_error('Mismatch in dimensions between potential and density', calledby='dftUPotential')

        IF(l_spinoffd.AND.spin_dim < 3) THEN
            CALL juDFT_error('Spin offdiagonal parts missing', calledby='dftUPotential')
        ENDIF

        IF(.NOT.l_spinoffd) spin_dim = MIN(2,spin_dim)


        CALL uj2f(jspins,ldau,f0,f2,f4,f6)

        ALLOCATE ( umatrix(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                           -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), source=0.0 )
        ALLOCATE(Vdc(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, SIZE(density,3)), source=cmplx_0)
        CALL umtx(ldau,f0,f2,f4,f6,umatrix)

        !--------------------------------------------------------!
        !  s       --                                        s'  !
        ! V     =  >  ( <m,p|V|m',q> - <m,p|V|q,m'> d     ) n    !
        !  m,m'    --                                s,s'    p,q !
        !        p,q,s'                                          !
        !--------------------------------------------------------!
        potential = cmplx_0

        DO ispin = 1, spin_dim
            DO m = -ldau%l,ldau%l
                DO mp =-ldau%l,ldau%l
                    vu = cmplx_0
                    IF(ispin < 3) THEN
                        DO jspin = 1, jspins
                        IF (ispin == jspin) THEN
                            DO p = -ldau%l,ldau%l
                                DO q = -ldau%l,ldau%l
                                    vu = vu + density(p, q,jspin) * ( umatrix(m,p,mp,q) - umatrix(m,p,q,mp) )
                                END DO
                            END DO
                        END IF
                        IF (ispin /= jspin .OR. jspins == 1) THEN
                            DO p = -ldau%l,ldau%l
                                DO q = -ldau%l,ldau%l
                                    vu = vu + umatrix(m,p,mp,q) * density(p, q,jspin)
                                END DO
                            END DO
                        END IF
                        END DO
                    ELSE
                        DO p = -ldau%l,ldau%l
                            DO q = -ldau%l,ldau%l
                                vu = vu - umatrix(m,p,q,mp) * density(p,q,ispin)
                            ENDDO
                        ENDDO
                    ENDIF
                    potential(m,mp,ispin) = potential(m,mp,ispin) + vu
                END DO ! m' loop
            END DO   ! m  loop
        END DO      ! outer spin loop

        !Take into account spin-degeneracy for non-spinpolarized calculations
        potential = potential * spin_deg

        !----------------------------------------------------------------------+
        !              s                                                       !
        !  ee      1  ---   s        s                     1        s  1       !
        ! E  (n) = -  >    n      ( V     + d     ( U (n - -) - J (n - -) ))   !
        !          2  ---   m,m'     m,m'    m,m'          2           2       !
        !             m,m'                                                     !
        !----------------------------------------------------------------------+

        interaction_energy = 0.0
        DO ispin = 1,spin_dim
            DO m = -ldau%l,ldau%l
                DO mp = -ldau%l,ldau%l
                    IF(ispin < 3) THEN
                        interaction_energy = interaction_energy + REAL(potential(m,mp,ispin)*density(m,mp,ispin))
                    ELSE
                        DO p = -ldau%l,ldau%l
                            DO q = -ldau%l,ldau%l
                                interaction_energy = interaction_energy + umatrix(m,p,q,mp) *&
                                                REAL( density(m,mp,ispin)*conjg(density(q,p,ispin)) &
                                                + conjg(density(mp,m,ispin))*density(p,q,ispin) )
                            ENDDO
                        ENDDO
                    ENDIF
                END DO
            END DO
        END DO
        interaction_energy = interaction_energy/2.0

    end subroutine

end module m_coulombPotential