MODULE m_crystalfield

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_crystalfield
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION:
   !>       -calculates the crystal-field-contrbution for the local hamiltonian
   !
   !------------------------------------------------------------------------------
   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_kkintgr
   USE m_ind_greensf

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE crystal_field(atoms,input,greensfCoeffs,hub1,vu)

      !calculates the crystal-field matrix for the local hamiltonian

      IMPLICIT NONE

      !-Type Arguments
      TYPE(t_greensfCoeffs), INTENT(IN)    :: greensfCoeffs
      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_hub1ham),       INTENT(INOUT) :: hub1

      !-Array Arguments
      COMPLEX,               INTENT(IN)    :: vu(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins) !LDA+U potential (should be removed from h_loc)

      !-Local Scalars
      INTEGER i_gf,l,nType,jspin,m,mp,ie,i_hia,kkcut,spin_cut

      !-Local Arrays
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      REAL :: integrand(greensfCoeffs%ne), norm(input%jspins),tr_hloc(atoms%n_hia,input%jspins)

      h_loc = 0.0
      tr_hloc = 0.0
      DO i_hia = 1, atoms%n_hia

         l     = atoms%lda_u(atoms%n_u+i_hia)%l
         nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
         i_gf = ind_greensf(atoms,l,nType)
         !Perform the integration 
         !
         ! \int_{E_b}^{E_c} dE E * N_LL'(E)
         !
         norm = 0.0
         DO jspin = 1, input%jspins
            spin_cut = MERGE(1,jspin,jspin.GT.2)
            kkcut = greensfCoeffs%kkintgr_cutoff(i_gf,spin_cut,2)
            DO m = -l, l
               DO mp = -l, l
                  integrand = 0.0
                  DO ie = 1, kkcut
                     integrand(ie) = ((ie-1) * greensfCoeffs%del+greensfCoeffs%e_bot) * greensfCoeffs%projdos(ie,i_gf,m,mp,jspin)
                  ENDDO
                  h_loc(m,mp,i_hia,jspin) = trapz(integrand,greensfCoeffs%del,greensfCoeffs%ne)
               ENDDO
               !trace of the integrated E*projdos
               tr_hloc(i_hia,jspin) = tr_hloc(i_hia,jspin) + h_loc(m,m,i_hia,jspin) - REAL(vu(m,m,i_hia,jspin))
            ENDDO
         ENDDO

         !Average over spins
         hub1%ccfmat(i_hia,:,:) = 0.0
         DO m = -l, l
            DO mp = -l, l
               hub1%ccfmat(i_hia,m,mp) = SUM(h_loc(m,mp,i_hia,:))/REAL(input%jspins) &
                                       - SUM(vu(m,mp,i_hia,:))/REAL(input%jspins)
            ENDDO
            hub1%ccfmat(i_hia,m,m) = hub1%ccfmat(i_hia,m,m) - SUM(tr_hloc(i_hia,:))/REAL(input%jspins*(2*l+1))
         ENDDO
      ENDDO


   END SUBROUTINE crystal_field

END MODULE m_crystalfield