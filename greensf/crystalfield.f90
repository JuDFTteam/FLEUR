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
   USE m_sgml
   USE m_anglso

   IMPLICIT NONE

   LOGICAL, PARAMETER :: l_debug = .FALSE.

   CONTAINS

   SUBROUTINE crystal_field(atoms,input,noco,greensfCoeffs,hub1,v)

      !calculates the crystal-field matrix for the local hamiltonian

      IMPLICIT NONE

      !-Type Arguments
      TYPE(t_greensfCoeffs), INTENT(IN)    :: greensfCoeffs
      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_noco),          INTENT(IN)    :: noco
      TYPE(t_hub1ham),       INTENT(INOUT) :: hub1

      !-Array Arguments
      TYPE(t_potden),        INTENT(IN)    :: v !LDA+U potential (should be removed from h_loc)

      !-Local Scalars
      INTEGER i_gf,l,nType,jspin,m,mp,ie,i_hia,kkcut,i_u,isp
      REAL    tr,xiSOC
      COMPLEX vso
      !-Local Arrays
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      REAL :: ex(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      REAL :: integrand(greensfCoeffs%ne), norm(greensfCoeffs%ne)

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
            kkcut = greensfCoeffs%kkintgr_cutoff(i_gf,jspin,2)
            norm = 0.0
            DO m = -l, l
               DO mp = -l, l
                  integrand = 0.0
                  DO ie = 1, kkcut
                     integrand(ie) = -1.0/pi_const * ((ie-1) * greensfCoeffs%del+greensfCoeffs%e_bot) &
                                     * REAL(greensfCoeffs%projdos(ie,m,mp,0,i_gf,jspin)/(3.0-input%jspins))
                     IF(m.EQ.mp) norm(ie) = norm(ie) -1.0/pi_const * REAL(greensfCoeffs%projdos(ie,m,mp,0,i_gf,jspin))/(3.0-input%jspins)
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
         !Remove LDA+U potential
         i_u = atoms%n_u+i_hia !position in the v%mmpmat array
         DO jspin = 1, input%jspins
            DO m = -l, l
               DO mp = -l, l
                  IF(ABS(h_loc(m,mp,i_hia,jspin)).GT.0.001.OR.m.EQ.mp) THEN
                     h_loc(m,mp,i_hia,jspin) = h_loc(m,mp,i_hia,jspin) - REAL(v%mmpmat(m,mp,i_u,jspin))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         !Remove SOC potential (only spin-diagonal)
         DO jspin = 1, 2
            DO m = -l, l
               DO mp = -l, l
                  isp = 3.0-2.0*jspin !1,-1
                  IF((ABS(noco%theta).LT.0.00001).AND.(ABS(noco%phi).LT.0.00001)) THEN
                     vso = CMPLX(sgml(l,m,isp,l,mp,isp),0.0)
                  ELSE
                     vso = anglso(noco%theta,noco%phi,l,m,isp,l,mp,isp)
                  ENDIF
                  h_loc(m,mp,i_hia,jspin) = h_loc(m,mp,i_hia,jspin) - vso/2.0 * hub1%xi(i_hia)/hartree_to_ev_const
               ENDDO
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "UP-REMOVED"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,1)
            WRITE(*,*) "DOWN-REMOVED"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,2)
         ENDIF
         ex = 0.0
         DO m= -l, l
            DO mp = -l, l
               ex(m,mp) = h_loc(m,mp,i_hia,1)-h_loc(m,mp,i_hia,2)
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "Exchange (eV)"
            WRITE(*,"(7f7.3)") ex(-3:3,-3:3)*hartree_to_ev_const*0.5
         ENDIF
         !------------------------------------------------------------------------------------
         ! If states move close to the cutoff we get some shift in the results on the diagonal
         ! The reason for this is a bit unclear we remove these results and replace them
         ! with either the -m or corresponding opposite spin result (only diagonal)
         !------------------------------------------------------------------------------------
        IF(.FALSE.) THEN
        DO m = -l, l
            !100 meV cutoff
            IF(ABS(ex(m,m)).LT.0.1/hartree_to_ev_const) CYCLE
            IF(ex(-m,-m).LT.0.1/hartree_to_ev_const) THEN
            !Assume the error is on the spin-down states
            IF(m.EQ.0) THEN
               WRITE(*,*) "Replacing m 0 spin down with m 0 spin up"
               h_loc(0,0,i_hia,2) = h_loc(0,0,i_hia,1)
            ELSE
               WRITE(*,*) "Replacing m ", m
               h_loc(m,m,i_hia,2) = h_loc(-m,-m,i_hia,2)
            ENDIF
         ENDIF
         ENDDO
        ENDIF

         !Average over spins
         hub1%ccfmat(i_hia,:,:) = 0.0
         DO m = -l, l
            DO mp = -l, l
               hub1%ccfmat(i_hia,m,mp) = SUM(h_loc(m,mp,i_hia,:))/2.0
               !For jspins.EQ.1 we need to take care of the fact that the spin-orbit coupling is opposite in spin 1/2
               IF(input%jspins.EQ.1) hub1%ccfmat(i_hia,m,mp) = (h_loc(m,mp,i_hia,1)+h_loc(-m,-mp,i_hia,1))/2.0
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "Average"
            WRITE(*,"(7f7.3)") hub1%ccfmat(i_hia,-3:3,-3:3)
         ENDIF
         DO m = -l, l
            DO mp = -l, l
               hub1%ccfmat(i_hia,m,mp) = (hub1%ccfmat(i_hia,m,mp)+hub1%ccfmat(i_hia,-m,-mp))/2.0
               hub1%ccfmat(i_hia,-m,-mp) = hub1%ccfmat(i_hia,m,mp)
            ENDDO
         ENDDO
         IF(l_debug) THEN
            WRITE(*,*) "SOC"
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

         IF(l_debug) THEN
            WRITE(*,*) "TRACELESS (eV)"
            WRITE(*,"(7f7.3)") hub1%ccfmat(i_hia,-3:3,-3:3)*hartree_to_ev_const
         ENDIF
      ENDDO
   END SUBROUTINE crystal_field

END MODULE m_crystalfield
