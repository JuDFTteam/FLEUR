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
   USE m_anglso

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE crystal_field(atoms,gfinp,hub1inp,input,nococonv,greensfImagPart,v,ef,hub1data)

      !calculates the crystal-field matrix for the local hamiltonian

      TYPE(t_greensfImagPart),   INTENT(IN)    :: greensfImagPart
      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_gfinp),             INTENT(IN)    :: gfinp
      TYPE(t_input),             INTENT(IN)    :: input
      TYPE(t_nococonv),          INTENT(IN)    :: nococonv
      TYPE(t_hub1inp),           INTENT(IN)    :: hub1inp
      TYPE(t_potden),            INTENT(IN)    :: v !LDA+U potential (should be removed from h_loc)
      REAL,                      INTENT(IN)    :: ef
      TYPE(t_hub1data),          INTENT(INOUT) :: hub1data

      !-Local Scalars
      INTEGER i_gf,l,nType,jspin,m,mp,ie,i_hia,kkcut,i_u,isp,iContour
      REAL    tr,xiSOC,del,eb
      COMPLEX vso
      !-Local Arrays
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      REAL :: ex(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      REAL :: integrand(gfinp%ne), norm(gfinp%ne)

      h_loc = 0.0
      DO i_hia = 1, atoms%n_hia

         l     = atoms%lda_u(atoms%n_u+i_hia)%l
         nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
         iContour = gfinp%hiaContour(i_hia)
         i_gf = gfinp%find(l,nType,iContour=iContour)
         !---------------------------------------------------------
         ! Perform the integration
         !---------------------------------------------------------
         ! \int_{E_b}^{E_c} dE E * N_LL'(E)
         !---------------------------------------------------------
         ! N_LL' is the l projected density of states and
         ! connected to the imaginary part of the greens function
         ! with a factor -1/pi
         !---------------------------------------------------------
         CALL gfinp%eMesh(ef,del,eb)
         DO jspin = 1, input%jspins
            !Use the same cutoffs as in the kramer kronigs integration
            kkcut = greensfImagPart%kkintgr_cutoff(i_gf,jspin,2)
            norm = 0.0
            DO m = -l, l
               DO mp = -l, l
                  integrand = 0.0
                  DO ie = 1, kkcut
                     integrand(ie) = -1.0/pi_const * ((ie-1) * del+eb) &
                                     * REAL(greensfImagPart%sphavg(ie,m,mp,i_gf,jspin)/(3.0-input%jspins))
                     IF(m.EQ.mp) norm(ie) = norm(ie) -1.0/pi_const * REAL(greensfImagPart%sphavg(ie,m,mp,i_gf,jspin))/(3.0-input%jspins)
                  ENDDO
                  h_loc(m,mp,i_hia,jspin) = trapz(integrand(1:kkcut),del,kkcut)
               ENDDO
            ENDDO
         ENDDO
#ifdef CPP_DEBUG
         WRITE(*,*) "UP"
         WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,1)
         IF(input%jspins.EQ.2) WRITE(*,*) "DOWN"
         IF(input%jspins.EQ.2) WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,2)
#endif
         !Remove LDA+U potential
         i_u = atoms%n_u+i_hia !position in the v%mmpmat array
         DO jspin = 1, input%jspins
            DO m = -l, l
               DO mp = -l, l
                  IF(ABS(REAL(v%mmpmat(m,mp,i_u,jspin))).GT.1e-4) h_loc(m,mp,i_hia,jspin) = h_loc(m,mp,i_hia,jspin) - REAL(v%mmpmat(m,mp,i_u,jspin))
               ENDDO
            ENDDO
         ENDDO
         !Remove SOC potential (only spin-diagonal)
         DO jspin = 1, 2
            DO m = -l, l
               DO mp = -l, l
                  isp = 3.0-2.0*jspin !1,-1
                  IF((ABS(nococonv%theta).LT.1e-5).AND.(ABS(nococonv%phi).LT.1e-5)) THEN
                     vso = CMPLX(sgml(l,m,isp,l,mp,isp),0.0)
                  ELSE
                     vso = anglso(nococonv%theta,nococonv%phi,l,m,isp,l,mp,isp)
                  ENDIF
                  h_loc(m,mp,i_hia,jspin) = h_loc(m,mp,i_hia,jspin) - REAL(vso)/2.0 * hub1data%xi(i_hia)/hartree_to_ev_const
               ENDDO
            ENDDO
         ENDDO
#ifdef CPP_DEBUG
         WRITE(*,*) "UP-REMOVED"
         WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,1)
         IF(input%jspins.EQ.2) WRITE(*,*) "DOWN-REMOVED"
         IF(input%jspins.EQ.2) WRITE(*,"(7f7.3)") h_loc(-3:3,-3:3,i_hia,2)
#endif
         ex = 0.0
         DO m= -l, l
            DO mp = -l, l
               ex(m,mp) = h_loc(m,mp,i_hia,1)-h_loc(m,mp,i_hia,2)
            ENDDO
         ENDDO
#ifdef CPP_DEBUG
         IF(input%jspins.EQ.2) WRITE(*,*) "Exchange (eV)"
         IF(input%jspins.EQ.2) WRITE(*,"(7f7.3)") ex(-3:3,-3:3)*hartree_to_ev_const*0.5
#endif
         !------------------------------------------------------------------------------------
         ! If states move close to the cutoff we get some shift in the results on the diagonal
         ! The reason for this is a bit unclear we remove these results and replace them
         ! with either the -m or corresponding opposite spin result (only diagonal)
         !------------------------------------------------------------------------------------
        !IF(.FALSE.) THEN
        !DO m = -l, l
            !100 meV cutoff
            !IF(ABS(ex(m,m)).LT.0.1/hartree_to_ev_const) CYCLE
            !IF(ex(-m,-m).LT.0.1/hartree_to_ev_const) THEN
            !Assume the error is on the spin-down states
            !IF(m.EQ.0) THEN
            !   WRITE(*,*) "Replacing m 0 spin down with m 0 spin up"
            !   h_loc(0,0,i_hia,2) = h_loc(0,0,i_hia,1)
            !ELSE
            !   WRITE(*,*) "Replacing m ", m
            !   h_loc(m,m,i_hia,2) = h_loc(-m,-m,i_hia,2)
            !ENDIF
         !ENDIF
         !ENDDO
        !ENDIF

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
         DO m = -l, l
            DO mp = -l, l
               hub1data%ccfmat(i_hia,m,mp) = (hub1data%ccfmat(i_hia,m,mp)+hub1data%ccfmat(i_hia,-m,-mp))/2.0
               hub1data%ccfmat(i_hia,-m,-mp) = hub1data%ccfmat(i_hia,m,mp)
            ENDDO
         ENDDO
#ifdef CPP_DEBUG
         WRITE(*,*) "SOC"
         WRITE(*,"(7f7.3)") hub1data%ccfmat(i_hia,-3:3,-3:3)
#endif
         tr = 0.0
         !calculate the trace
         DO m = -l, l
            tr = tr + hub1data%ccfmat(i_hia,m,m)
         ENDDO
#ifdef CPP_DEPUG
         WRITE(*,*) "TRACE"
         WRITE(*,"(2f7.3)") tr, tr/(2*l+1)
#endif
         !Remove trace
         DO m = -l, l
            hub1data%ccfmat(i_hia,m,m) = hub1data%ccfmat(i_hia,m,m) - tr/(2*l+1)
         ENDDO

#ifdef CPP_DEPUG
         WRITE(*,*) "TRACELESS (eV)"
         WRITE(*,"(7f7.3)") hub1data%ccfmat(i_hia,-3:3,-3:3)*hartree_to_ev_const
#endif
      ENDDO
   END SUBROUTINE crystal_field

END MODULE m_crystalfield
