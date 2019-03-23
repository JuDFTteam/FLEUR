MODULE m_eff_excinteraction

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_eff_excinteraction
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  This module calculates the effective exchange interaction from th onsite
   !>  green's function according to Condens. Matter 26 (2014) 476003 EQ.1
   !
   ! REVISION HISTORY:
   ! February 2019 - Initial Version
   ! March    2019 - Changed calculation of the onsite exchange matrix
   !------------------------------------------------------------------------------

   CONTAINS

   SUBROUTINE eff_excinteraction(gOnsite,atoms,input,j0,onsite_exc_split)

      USE m_types
      USE m_constants
      USE m_juDFT
      

      IMPLICIT NONE

      TYPE(t_greensf),        INTENT(IN)  :: gOnsite
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      REAL,                   INTENT(OUT) :: j0
      TYPE(t_input),          INTENT(IN)  :: input
      REAL,                   INTENT(IN)  :: onsite_exc_split

      COMPLEX tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),integrand(2)
      INTEGER i_j0,iz,m,l,mp,ispin,n,i_gf,matsize,i
      LOGICAL l_matinv



      COMPLEX :: delta(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX :: g_up(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX :: g_dwn(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      INTEGER              :: info
      INTEGER, ALLOCATABLE :: ipiv(:)
      COMPLEX, ALLOCATABLE :: work(:)
      COMPLEX, ALLOCATABLE :: invup(:,:)
      COMPLEX, ALLOCATABLE :: invdwn(:,:)
      COMPLEX, ALLOCATABLE :: inv(:,:)

      l_matinv = .true. !Determines how the onsite exchange splitting is calculated
      
      DO i_j0 = 1, atoms%n_j0
         j0 = 0.0
         l = atoms%j0(i_j0)%l
         n = atoms%j0(i_j0)%atomType
         !Find the corresponding index of the onsite gf
         CALL gOnsite%index(l,n,i_gf)
         IF(.true.) THEN
            matsize = 2*l+1
            ALLOCATE(work(matsize))
            ALLOCATE(invup(matsize,matsize))
            ALLOCATE(invdwn(matsize,matsize))
            ALLOCATE(inv(matsize,matsize))
            ALLOCATE(ipiv(matsize))
         ENDIF 

         DO iz = 1, gOnsite%nz
            !
            !calculate the onsite exchange matrix
            !
            IF(l_matinv) THEN
               !First Way: Matrix Inversion
               !---------------------------------------------
               !\Delta = (G_up)^-1-(G_down)^-1
               !---------------------------------------------
               !Symmetrize the green's function for up/down 
               !spin with respect to the complex plane
               !Here we assume that the onsite Hamiltonian
               !is real
               !---------------------------------------------
               !G^(up/down)^-1 = 1/2 * (G+^(up/down) + G-^(up/down))
               !---------------------------------------------
               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,1,1)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invup(1:matsize,1:matsize) = inv(1:matsize,1:matsize)
               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,1,2)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invup(1:matsize,1:matsize) = 1/2.0*(invup(1:matsize,1:matsize) +inv(1:matsize,1:matsize))


               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,2,1)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invdwn(1:matsize,1:matsize) = inv(1:matsize,1:matsize)
               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,2,2)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invdwn(1:matsize,1:matsize) = 1/2.0*(invdwn(1:matsize,1:matsize) +inv(1:matsize,1:matsize))

               delta(-l:l,-l:l) = invup(1:matsize,1:matsize) - invdwn(1:matsize,1:matsize)
            ELSE
               !Second Way: onsite_exc_split is the difference in the center of gravity of the up/down bands
               delta = 0.0
               DO m = -l, l
                  delta(m,m) = onsite_exc_split
               ENDDO
            ENDIF
            !
            !  Tr[\Delta (G_up-G-down) + \Delta G_up \Delta G-down]
            !
            integrand = 0.0
            DO i = 1, 2
               g_up(-l:l,-l:l)   = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,1,i)
               g_dwn(-l:l,-l:l)  = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,2,i)
               tmp(-l:l,-l:l)    = g_up(-l:l,-l:l)-g_dwn(-l:l,-l:l)
               tmp(-l:l,-l:l)    = matmul(delta(-l:l,-l:l),tmp(-l:l,-l:l))
               tmp(-l:l,-l:l)    = tmp(-l:l,-l:l) + matmul(matmul(delta(-l:l,-l:l),g_up(-l:l,-l:l)),&
                                 matmul(delta(-l:l,-l:l),g_dwn(-l:l,-l:l)))
               !Trace over m
               DO m = -l,l
                  integrand(i) = integrand(i) + tmp(m,m)
               ENDDO
            ENDDO
            j0 = j0 + 1/2.0 * AIMAG(integrand(1)*gOnsite%de(iz)-integrand(2)*conjg(gOnsite%de(iz)))
         ENDDO

         
         j0 = -j0*1/fpi_const*hartree_to_ev_const
         WRITE(*,*)  "Eff. Exchange Interaction for atom", n, ": ", j0, "eV"
         IF(ALLOCATED(work)) DEALLOCATE(work)
         IF(ALLOCATED(ipiv)) DEALLOCATE(ipiv)
         IF(ALLOCATED(invup)) DEALLOCATE(invup)
         IF(ALLOCATED(invdwn)) DEALLOCATE(invdwn)
      ENDDO
   END SUBROUTINE eff_excinteraction


END MODULE m_eff_excinteraction