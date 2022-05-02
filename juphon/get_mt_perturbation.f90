!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_get_mt_perturbation
    USE m_juDFT
    USE m_polangle
    USE m_types
    USE m_constants
    USE m_mt_tofrom_grid

    IMPLICIT NONE

CONTAINS

    SUBROUTINE get_mt_local_perturbation(atoms,sphhar,sym,noco,den,den1,den1im)
        TYPE(t_atoms),  INTENT(IN)    :: atoms
        TYPE(t_sphhar), INTENT(IN)    :: sphhar
        TYPE(t_sym),    INTENT(IN)    :: sym
        TYPE(t_noco),   INTENT(IN)    :: noco
        TYPE(t_potden), INTENT(INOUT) :: den, den1, den1im

        TYPE(t_gradients)             :: grad

        INTEGER                       :: n, nsp, imesh, i
        REAL                          :: rho_11, rho_22, m
        REAL                          :: rhotot, rho_up, rho_down, theta, phi
        COMPLEX                       :: n1, mx1, my1, mz1, m1, t1, p1, rho1_up, rho1_down
        REAL                          :: eps=1E-10
        REAL, ALLOCATABLE             :: ch(:,:), chre(:,:), chim(:,:)

        nsp=atoms%nsp()
        ALLOCATE(ch(nsp*atoms%jmtd,4), den1%theta_mt(nsp*atoms%jmtd,atoms%ntype), &
                                       den1%phi_mt(nsp*atoms%jmtd,atoms%ntype))
        ALLOCATE(chre(nsp*atoms%jmtd,4), chim(nsp*atoms%jmtd,4))
        ALLOCATE(den1im%theta_mt(nsp*atoms%jmtd,atoms%ntype), den1im%phi_mt(nsp*atoms%jmtd,atoms%ntype))

        CALL init_mt_grid(4,atoms,sphhar,.FALSE.,sym)

        DO n=1,atoms%ntype
            CALL mt_to_grid(.FALSE., 2, atoms, sym, sphhar, .FALSE., den%mt(:,0:,n,:), n, noco, grad, ch)
            CALL mt_to_grid(.FALSE., 4, atoms, sym, sphhar, .FALSE., den1%mt(:,0:,n,:), n, noco, grad, chre)
            CALL mt_to_grid(.FALSE., 4, atoms, sym, sphhar, .FALSE., den1im%mt(:,0:,n,:), n, noco, grad, chim)

            DO imesh = 1, nsp*atoms%jri(n)
                rho_11   = ch(imesh,1)
                rho_22   = ch(imesh,2)

                m = rho_11 - rho_22

                ! Calculate perturbed total and magnetization density
                n1  = chre(imesh,1) + ImagUnit * chim(imesh,1)
                mx1 = chre(imesh,2) + ImagUnit * chim(imesh,2)
                my1 = chre(imesh,3) + ImagUnit * chim(imesh,3)
                mz1 = chre(imesh,4) + ImagUnit * chim(imesh,4)

                theta = den%theta_mt(imesh,n)
                phi   = den%phi_mt(imesh,n)

                ! Calculate the perturbed absolute value of the magnetization density
                m1 = cmplx(0.0, 0.0)
                m1 = m1 + mx1 * SIN(theta) * COS(phi)
                m1 = m1 + my1 * SIN(theta) * SIN(phi)
                m1 = m1 + mz1 * COS(theta)

                ! Calculate the perturbed angles
                t1 = cmplx(0.0, 0.0)
                t1 = t1 + mx1 * COS(theta) * COS(phi) / m
                t1 = t1 + my1 * COS(theta) * SIN(phi) / m
                t1 = t1 - mz1 * SIN(theta)            / m

                p1 = cmplx(0.0, 0.0)
                p1 = p1 - mx1 * SIN(theta) * SIN(phi) / m
                p1 = p1 + my1 * SIN(theta) * COS(phi) / m

                rho1_up   = (n1 + m1)/2
                rho1_down = (n1 - m1)/2

                chre(imesh,1) =  REAL(rho1_up  )
                chre(imesh,2) =  REAL(rho1_down)
                chim(imesh,1) = AIMAG(rho1_up  )
                chim(imesh,2) = AIMAG(rho1_down)
                den1%theta_mt(imesh,n)   =  REAL(t1)
                den1%phi_mt(imesh,n)     =  REAL(p1)
                den1im%theta_mt(imesh,n) = AIMAG(t1)
                den1im%phi_mt(imesh,n)   = AIMAG(p1)
            END DO
            den1%mt(:,0:,n,:)=0.0
            den1im%mt(:,0:,n,:)=0.0
            CALL mt_from_grid(atoms,sym,sphhar,n,2,chre,den1%mt(:,0:,n,:))
            CALL mt_from_grid(atoms,sym,sphhar,n,2,chim,den1im%mt(:,0:,n,:))
            DO i=1,atoms%jri(n)
                den1%mt(i,:,n,:)=den1%mt(i,:,n,:)*atoms%rmsh(i,n)**2
                den1im%mt(i,:,n,:)=den1im%mt(i,:,n,:)*atoms%rmsh(i,n)**2
            END DO
        END DO

        CALL finish_mt_grid()

    END SUBROUTINE get_mt_local_perturbation

    SUBROUTINE get_mt_global_perturbation(atoms,sphhar,sym,den,den1,den1im,noco,vtot,vtot1,vtot1im)
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_potden), INTENT(IN)    :: den, den1, den1im
      TYPE(t_noco),   INTENT(IN)    :: noco
      TYPE(t_potden), INTENT(INOUT) :: vtot,vtot1,vtot1im

      TYPE(t_gradients)             :: grad

      INTEGER                       :: n, nsp, imesh, i
      REAL                          :: theta, phi, v11, v22
      COMPLEX                       :: v1up, v1down, v1eff, b1eff
      REAL, ALLOCATABLE             :: ch(:,:), chtmp(:,:), chretmp(:,:), chimtmp(:,:)

      COMPLEX :: a11, a22, a21, a12, av11, av22, av21, av12
      COMPLEX :: t1, p1, v21, v12, v1mat11, v1mat22, v1mat21, v1mat12

      nsp=atoms%nsp()
      ALLOCATE(chtmp(nsp*atoms%jmtd,4))
      ALLOCATE(chretmp(nsp*atoms%jmtd,4))
      ALLOCATE(chimtmp(nsp*atoms%jmtd,4))
      !TODO: Make sure the indices for rho1 are 1,2,3,4 == n1,mx1,my1,mz1
      CALL init_mt_grid(4,atoms,sphhar,.FALSE.,sym)
      DO n=1,atoms%ntype

         DO i=1,atoms%jri(n)
            vtot%mt(i,:,n,:)=vtot%mt(i,:,n,:)*atoms%rmsh(i,n)**2
            vtot1%mt(i,:,n,:)=vtot1%mt(i,:,n,:)*atoms%rmsh(i,n)**2
            vtot1im%mt(i,:,n,:)=vtot1im%mt(i,:,n,:)*atoms%rmsh(i,n)**2
         END DO

         CALL mt_to_grid(.FALSE.,4,atoms,sym,sphhar,.FALSE.,vtot%mt(:,0:,n,:),n,noco,grad,chtmp(:,1:2))
         CALL mt_to_grid(.FALSE.,2,atoms,sym,sphhar,.FALSE.,vtot1%mt(:,0:,n,:),n,noco,grad,chretmp(:,1:2))
         CALL mt_to_grid(.FALSE.,2,atoms,sym,sphhar,.FALSE.,vtot1im%mt(:,0:,n,:),n,noco,grad,chimtmp(:,1:2))

         DO imesh = 1, nsp*atoms%jri(n)
             v11   = chtmp(imesh, 1)
             v22   = chtmp(imesh, 2)
             v21   = chtmp(imesh, 3) + ImagUnit * chtmp(imesh, 4)
             v12   = chtmp(imesh, 3) - ImagUnit * chtmp(imesh, 4)

            v1up         = chretmp(imesh,1) + ImagUnit * chimtmp(imesh,1)
            v1down       = chretmp(imesh,2) + ImagUnit * chimtmp(imesh,2)

            theta       = den%theta_mt(imesh,n)
            phi         = den%phi_mt(imesh,n)

            v1eff        = (v1up + v1down)/2.0
            b1eff        = (v1up - v1down)/2.0

            v1mat11 = v1eff + b1eff * COS(theta)                    !11
            v1mat22 = v1eff - b1eff * COS(theta)                    !22
            v1mat21 =         b1eff * SIN(theta)*EXP( Imagunit*phi) !21
            v1mat12 =         b1eff * SIN(theta)*EXP(-Imagunit*phi) !12

            t1 = den1%theta_mt(imesh,n) + ImagUnit * den1im%theta_mt(imesh,n)
            p1 = den1%phi_mt(imesh,n) + ImagUnit * den1im%phi_mt(imesh,n)

            a11 =      -ImagUnit      * p1 / 2.0
            a22 =       ImagUnit      * p1 / 2.0
            a21 =  EXP( ImagUnit*phi) * t1 / 2.0
            a12 = -EXP(-ImagUnit*phi) * t1 / 2.0

            av11 = a11 * v11 + a12 * v21 !11
            av22 = a21 * v12 + a22 * v22 !22
            av21 = a21 * v11 + a22 * v21 !21
            av12 = a11 * v12 + a12 * v22 !12

            v1mat11 = v1mat11 + av11 + CONJG(av11)
            v1mat22 = v1mat22 + av22 + CONJG(av22)
            v1mat21 = v1mat21 + av21 + CONJG(av12)
            v1mat12 = v1mat12 + av12 + CONJG(av21)

            chretmp(imesh, 1) =  REAL(v1mat11)
            chretmp(imesh, 2) =  REAL(v1mat22)
            chretmp(imesh, 3) =  REAL(v1mat21)
            chretmp(imesh, 4) =  REAL(v1mat12)

            chimtmp(imesh, 1) = AIMAG(v1mat11)
            chimtmp(imesh, 2) = AIMAG(v1mat22)
            chimtmp(imesh, 3) = AIMAG(v1mat21)
            chimtmp(imesh, 4) = AIMAG(v1mat12)

         END DO

         vtot1%mt(:,0:,n,:)=0.0
         vtot1im%mt(:,0:,n,:)=0.0

         CALL mt_from_grid(atoms,sym,sphhar,n,4,chretmp,vtot1%mt(:,0:,n,:))
         CALL mt_from_grid(atoms,sym,sphhar,n,4,chimtmp,vtot1im%mt(:,0:,n,:))

      END DO

      CALL finish_mt_grid()

    END SUBROUTINE get_mt_global_perturbation

END MODULE m_get_mt_perturbation
