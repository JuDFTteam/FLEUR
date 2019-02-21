MODULE m_eff_excinteraction

   !This module calculates the effective exchange interaction
   !from the onsite green's function

   CONTAINS

   SUBROUTINE eff_excinteraction(gOnsite,atoms,input,j0)

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
      INTEGER i_j0,iz,m,l,mp,ispin,n,i_gf,ind

      REAL f0(atoms%n_j0,input%jspins),f2(atoms%n_j0,input%jspins)
      REAL f4(atoms%n_j0,input%jspins),f6(atoms%n_j0,input%jspins)
      REAL, ALLOCATABLE :: u(:,:,:,:,:)
      REAL, ALLOCATABLE :: delta(:,:)

      ! calculate slater integrals from u and j
      CALL uj2f(input%jspins,atoms%j0(:),atoms%n_j0,f0,f2,f4,f6)

      ! set up e-e- interaction matrix
      ALLOCATE (u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_j0)))
      ALLOCATE (delta(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))

      !If the slater integrals are calculated in cored they depend on spin
      !We average over them
      f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
      f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
      f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
      f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2

      CALL umtx(atoms%j0(:),atoms%n_j0,f0(:,1),f2(:,1),f4(:,1),f6(:,1),u)

      DO i_j0 = 1, atoms%n_j0
         !
         !Calculate the onsite exchange matrix 
         !
         j0 = 0.0
         l = atoms%j0(i_j0)%l
         n = atoms%j0(i_j0)%atomType
         delta = 0.0
         DO m = -l,l
            DO mp = -l,l
               IF(m.NE.mp) delta(m,mp) = u(m,mp,mp,m,i_j0)
            ENDDO
         ENDDO

         WRITE(*,'(7f14.8)') delta(:,:)

         !Find the corresponding index of the onsite gf
         ind =  0
         DO i_gf = 1, gOnsite%n_gf
            IF(gOnsite%atomType(i_gf).EQ.n.AND.gOnsite%l_gf(i_gf).EQ.l) THEN
               ind = i_gf
               EXIT
            ENDIF
         ENDDO

         DO iz = 1, gOnsite%nz
            !
            !calculate the integrand for the current energy point
            !
            !  Tr[\Delta (G_up-G-down) + \Delta G_up \Delta G-down]
            !
            tmp(-l:l,-l:l) = gOnsite%gmmpMat(1,iz,ind,-l:l,-l:l,2)-gOnsite%gmmpMat(1,iz,ind,-l:l,-l:l,1)
            tmp(-l:l,-l:l) = matmul(delta(-l:l,-l:l),tmp(-l:l,-l:l))
            tmp(-l:l,-l:l) = tmp(-l:l,-l:l) + matmul(matmul(delta(-l:l,-l:l),gOnsite%gmmpMat(1,iz,ind,-l:l,-l:l,2)),&
                              matmul(delta(-l:l,-l:l),gOnsite%gmmpMat(1,iz,ind,-l:l,-l:l,1)))
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

         WRITE(*,*)  "Eff. Exchange Interaction for atom", n, ": ", j0, "eV"
      ENDDO
   END SUBROUTINE eff_excinteraction

END MODULE m_eff_excinteraction