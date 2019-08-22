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

   LOGICAL, PARAMETER :: l_debug = .TRUE.

   CONTAINS

   SUBROUTINE crystal_field(atoms,input,greensfCoeffs,hub1,v)

      !calculates the crystal-field matrix for the local hamiltonian

      IMPLICIT NONE

      !-Type Arguments
      TYPE(t_greensfCoeffs), INTENT(IN)    :: greensfCoeffs
      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_hub1ham),       INTENT(INOUT) :: hub1

      !-Array Arguments
      TYPE(t_potden),        INTENT(IN)    :: v !LDA+U potential (should be removed from h_loc)

      !-Local Scalars
      INTEGER i_gf,l,nType,jspin,m,mp,ie,i_hia,kkcut,spin_cut,i_u
      REAL    tr,xiSOC
      !-Local Arrays
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      REAL :: integrand(greensfCoeffs%ne)

      h_loc = 0.0
      DO i_hia = 1, atoms%n_hia

         l     = atoms%lda_u(atoms%n_u+i_hia)%l
         nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
         i_gf = ind_greensf(atoms,l,nType)
         !---------------------------------------------------------
         ! Perform the integration 
         !---------------------------------------------------------
         ! \int_{E_b}^{E_c} dE E * N_LL'(E)
         !---------------------------------------------------------
         ! N_LL' is the l projected density of states and 
         ! connected to the imaginary part of the greens function 
         ! with a factor -1/pi
         !---------------------------------------------------------
         DO jspin = 1, input%jspins
            !Use the same cutoffs as in the kramer kronigs integration
            spin_cut = MERGE(1,jspin,jspin.GT.2)
            kkcut = greensfCoeffs%kkintgr_cutoff(i_gf,spin_cut,2)
            DO m = -l, l
               DO mp = -l, l
                  integrand = 0.0
                  DO ie = 1, kkcut
                     integrand(ie) = -1.0/pi_const * ((ie-1) * greensfCoeffs%del+greensfCoeffs%e_bot) &
                                     * greensfCoeffs%projdos(ie,i_gf,m,mp,jspin)/(3.0-input%jspins)
                  ENDDO
                  h_loc(m,mp,i_hia,jspin) = trapz(integrand(1:kkcut),greensfCoeffs%del,kkcut)
               ENDDO
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "UP"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,1)
            WRITE(*,*) "DOWN"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,2)
         ENDIF
         !Remove LDA+U and SOC potential 
         i_u = atoms%n_u+i_hia !position in the v%mmpmat array 
         DO jspin = 1, input%jspins
            DO m = -l, l
               DO mp = -l, l
                  !LDA+U potential
                  h_loc(m,mp,i_hia,jspin) = h_loc(m,mp,i_hia,jspin) - REAL(v%mmpmat(m,mp,i_u,jspin))
               ENDDO
               !h_loc(m,m,i_hia,jspin) = h_loc(m,m,i_hia,jspin) - hub1%xi(i_hia)/hartree_to_ev_const * m * (1.5-jspin) * MERGE(-1,1,input%jspins.EQ.1)
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "UP-REMOVED"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,1)
            WRITE(*,*) "DOWN-REMOVED"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,2)
         ENDIF
         !Average over spins
         hub1%ccfmat(i_hia,:,:) = 0.0
         DO m = -l, l
            DO mp = -l, l
               hub1%ccfmat(i_hia,m,mp) = SUM(h_loc(m,mp,i_hia,:))/2.0
               IF(input%jspins.EQ.1) hub1%ccfmat(i_hia,m,mp) = (h_loc(m,mp,i_hia,1)+h_loc(-m,-mp,i_hia,1))/2.0
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "Average"
            WRITE(*,"(7f7.3)") hub1%ccfmat(i_hia,-3:3,-3:3)
         ENDIF
         tr = 0.0
         !calculate the trace
         DO m = -l, l 
            tr = tr + hub1%ccfmat(i_hia,m,m)
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "TRACE"
            WRITE(*,"(2f7.3)") tr, tr/(2*l+1)
         ENDIF
         !Remove trace 
         DO m = -l, l 
            hub1%ccfmat(i_hia,m,m) = hub1%ccfmat(i_hia,m,m) - tr/(2*l+1) 
         ENDDO
         DO m = -l, l
            DO mp = -l, l
               hub1%ccfmat(i_hia,m,mp) = (hub1%ccfmat(i_hia,m,mp)+hub1%ccfmat(i_hia,-m,-mp))/2.0
               hub1%ccfmat(i_hia,-m,-mp) = hub1%ccfmat(i_hia,m,mp)
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "TRACELESS (eV)"
            WRITE(*,"(7f7.3)") hub1%ccfmat(i_hia,-3:3,-3:3)*hartree_to_ev_const
         ENDIF
      ENDDO
   END SUBROUTINE crystal_field

END MODULE m_crystalfield