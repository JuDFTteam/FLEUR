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
   USE m_trapz
   USE m_sgml
   USE m_rotMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE crystal_field(atoms,gfinp,input,noco,nococonv,greensFunction,v,ef,hub1data)

      !calculates the crystal-field matrix for the local hamiltonian

      TYPE(t_greensf),           INTENT(IN)    :: greensFunction(:)
      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_gfinp),             INTENT(IN)    :: gfinp
      TYPE(t_input),             INTENT(IN)    :: input
      TYPE(t_noco),              INTENT(IN)    :: noco
      TYPE(t_nococonv),          INTENT(IN)    :: nococonv
      TYPE(t_potden),            INTENT(IN)    :: v !LDA+U potential (should be removed from h_loc)
      REAL,                      INTENT(IN)    :: ef
      TYPE(t_hub1data),          INTENT(INOUT) :: hub1data

      INTEGER i_gf,l,atomType,jspin,m,mp,ie,i_hia,i_u,isp
      REAL    tr,del,eb
      COMPLEX vso
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX, ALLOCATABLE :: potmmpmat(:,:,:)


      ALLOCATE(potmmpmat(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, input%jspins))

      h_loc = 0.0
      DO i_hia = 1, atoms%n_hia

         if (gfinp%hiaFitElem(i_hia) < 0) cycle

         l     = atoms%lda_u(atoms%n_u+i_hia)%l
         atomType = atoms%lda_u(atoms%n_u+i_hia)%atomType

         i_gf = gfinp%hiaFitElem(i_hia)

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
            h_loc(:,:,i_hia,jspin) = greensFunction(i_gf)%spectralDensityMoment(1, jspin, gfinp, input, atoms, noco, nococonv)
         ENDDO

#ifdef CPP_DEBUG
         WRITE(*,*) "UP"
         WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,1)
         IF(input%jspins.EQ.2) THEN
            WRITE(*,*) "DOWN"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,2)
         ENDIF
#endif
         !Remove LDA+U potential
         i_u = atoms%n_u+i_hia !position in the v%mmpmat array
         !Rotate the occupation matrix into the global frame in real-space
         IF(noco%l_noco) THEN
            potmmpmat = rotMMPmat(v%mmpmat(:,:,i_u,:input%jspins),nococonv%alph(atomType),nococonv%beta(atomType),0.0,l, inverse=.true.)
         ELSE IF(noco%l_soc) THEN
            potmmpmat = rotMMPmat(v%mmpmat(:,:,i_u,:input%jspins),nococonv%phi,nococonv%theta,0.0,l, inverse=.true.)
         ENDIF
         DO jspin = 1, input%jspins
            DO m = -l, l
               DO mp = -l, l
                  IF(ABS(potmmpmat(m,mp,jspin)).GT.1e-4) THEN
                     h_loc(m,mp,i_hia,jspin) = h_loc(m,mp,i_hia,jspin) - potmmpmat(m,mp,jspin)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         IF(ABS(hub1data%xi(i_hia)).GT.1e-12) THEN
            !Remove SOC potential (only spin-diagonal)
            DO jspin = 1, 2
               DO m = -l, l
                  DO mp = -l, l
                     isp = 3-2*jspin !1,-1
                     vso = CMPLX(sgml(l,m,isp,l,mp,isp),0.0)
                     h_loc(m,mp,i_hia,jspin) = h_loc(m,mp,i_hia,jspin) - REAL(vso)/2.0 * hub1data%xi(i_hia)/hartree_to_ev_const
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
#ifdef CPP_DEBUG
         WRITE(*,*) "UP-REMOVED"
         WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,1)
         IF(input%jspins.EQ.2) THEN
            WRITE(*,*) "DOWN-REMOVED"
            WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,2)
         ENDIF
#endif

         !Average over spins
         hub1data%ccfmat(i_hia,:,:) = 0.0
         DO m = -l, l
            DO mp = -l, l
               hub1data%ccfmat(i_hia,m,mp) = SUM(h_loc(m,mp,i_hia,:))/2.0
               !For jspins.EQ.1 we need to take care of the fact that the spin-orbit coupling is opposite in spin 1/2
               IF(input%jspins.EQ.1) hub1data%ccfmat(i_hia,m,mp) = (h_loc(m,mp,i_hia,1)+h_loc(-m,-mp,i_hia,1))/2.0
            ENDDO
         ENDDO
#ifdef CPP_DEBUG
         WRITE(*,*) "Average"
         WRITE(*,"(7f7.3)") hub1data%ccfmat(i_hia,-3:3,-3:3)
#endif

         !-----------------------------------
         ! Symmetrize Matrix
         ! Delta_CF(m,mp) = Delta_CF(-m,-mp)
         !-----------------------------------
         IF(input%jspins.EQ.2) THEN
            DO m = -l, l
               DO mp = -l, l
                  hub1data%ccfmat(i_hia,m,mp) = (hub1data%ccfmat(i_hia,m,mp)+hub1data%ccfmat(i_hia,-m,-mp))/2.0
                  hub1data%ccfmat(i_hia,-m,-mp) = hub1data%ccfmat(i_hia,m,mp)
               ENDDO
            ENDDO
         ENDIF
#ifdef CPP_DEBUG
         WRITE(*,*) "Average symmetrized"
         WRITE(*,"(7f7.3)") hub1data%ccfmat(i_hia,-3:3,-3:3)
#endif

         tr = 0.0
         !calculate the trace
         DO m = -l, l
            tr = tr + hub1data%ccfmat(i_hia,m,m)
         ENDDO
#ifdef CPP_DEBUG
         WRITE(*,*) "TRACE"
         WRITE(*,"(2f7.3)") tr, tr/(2*l+1)
#endif
         !Remove trace
         DO m = -l, l
            hub1data%ccfmat(i_hia,m,m) = hub1data%ccfmat(i_hia,m,m) - tr/(2*l+1)
         ENDDO

#ifdef CPP_DEBUG
         WRITE(*,*) "TRACELESS ccf (eV)"
         WRITE(*,"(7f7.3)") hub1data%ccfmat(i_hia,-3:3,-3:3)*hartree_to_ev_const
#endif
      ENDDO
   end subroutine crystal_field

END MODULE m_crystalfield
