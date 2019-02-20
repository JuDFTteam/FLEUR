MODULE m_eff_excsplitting

   !This module calculates the effective exchange interaction
   !from the onsite green's function

   CONTAINS

   SUBROUTINE eff_excsplitting(gOnsite,atoms,input,j0)

      USE m_types
      USE m_constants
      USE m_umtx
      USE m_uj2f

      IMPLICIT NONE

      TYPE(t_greensf),        INTENT(IN)  :: gOnsite
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      REAL,                   INTENT(OUT) :: j0
      TYPE(t_input),          INTENT(IN)  :: input

      COMPLEX tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),integrand
      REAL f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins),f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      REAL, ALLOCATABLE :: u(:,:,:,:,:,:)
      REAL, ALLOCATABLE :: delta(:,:)

      INTEGER i_hia,iz,m,l,mp,ispin
      ! calculate slater integrals from u and j
      CALL uj2f(input%jspins,atoms%lda_hia(:),atoms%n_hia,f0,f2,f4,f6)

      ! set up e-e- interaction matrix
      ALLOCATE (u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_hia),input%jspins))
      ALLOCATE (delta(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
      DO ispin = 1, 1 ! input%jspins
         f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
         f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
         f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
         f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
         CALL umtx(atoms%lda_hia(:),atoms%n_hia,f0(1,ispin),f2(1,ispin),f4(1,ispin),f6(1,ispin),&
                 u(-lmaxU_const,-lmaxU_const,-lmaxU_const,-lmaxU_const,1,ispin))
      END DO

      DO i_hia = 1, atoms%n_hia
         !
         !Calculate the onsite exchange matrix 
         !
         j0 = 0.0
         l = atoms%lda_hia(i_hia)%l
         delta = 0.0
         DO m = -l,l
            DO mp = -l,l
               IF(m.NE.mp) delta(m,mp) = u(m,mp,mp,m,atoms%n_u+i_hia,1)
            ENDDO
         ENDDO

         WRITE(*,'(7f14.8)') delta(:,:)

         DO iz = 1, gOnsite%nz
            !
            !calculate the integrand for the current energy point
            !
            !  Tr[\Delta (G_up-G-down) + \Delta G_up \Delta G-down]
            !
            tmp(-l:l,-l:l) = gOnsite%gmmpMat(1,iz,i_hia,-l:l,-l:l,2)-gOnsite%gmmpMat(1,iz,i_hia,-l:l,-l:l,1)
            tmp(-l:l,-l:l) = matmul(delta(-l:l,-l:l),tmp(-l:l,-l:l))
            tmp(-l:l,-l:l) = tmp(-l:l,-l:l) + matmul(matmul(delta(-l:l,-l:l),gOnsite%gmmpMat(1,iz,i_hia,-l:l,-l:l,2)),&
                              matmul(delta(-l:l,-l:l),gOnsite%gmmpMat(1,iz,i_hia,-l:l,-l:l,1)))
            !WRITE(*,*) iz
            !WRITE(*,'(14f14.8)') tmp(:,:)
            integrand = CMPLX(0.0,0.0)
            !Trace over m
            DO m = -l,l
               integrand = integrand + tmp(m,m)
            ENDDO

            j0 = j0 + AIMAG(integrand*gOnsite%de(iz))
         ENDDO

         j0 = j0*1/fpi_const*hartree_to_ev_const

         WRITE(*,*)  "Eff. Exchange splitting for atom", atoms%lda_hia(i_hia)%atomType, ": ", j0, "eV"
      ENDDO
   END SUBROUTINE eff_excsplitting

END MODULE m_eff_excsplitting