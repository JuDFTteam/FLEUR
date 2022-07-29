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

   SUBROUTINE crystal_field(atoms,gfinp,input,noco,nococonv,greensfImagPart,v,ef,hub1data)

      !calculates the crystal-field matrix for the local hamiltonian

      TYPE(t_greensfImagPart),   INTENT(IN)    :: greensfImagPart
      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_gfinp),             INTENT(IN)    :: gfinp
      TYPE(t_input),             INTENT(IN)    :: input
      TYPE(t_noco),              INTENT(IN)    :: noco
      TYPE(t_nococonv),          INTENT(IN)    :: nococonv
      TYPE(t_potden),            INTENT(IN)    :: v !LDA+U potential (should be removed from h_loc)
      REAL,                      INTENT(IN)    :: ef
      TYPE(t_hub1data),          INTENT(INOUT) :: hub1data

      !-Local Scalars
      INTEGER i_gf,l,nType,jspin,m,mp,ie,i_hia,i_u,isp,i_elem
      REAL    tr,del,eb
      COMPLEX vso
      LOGICAL, PARAMETER :: l_correctMinus = .FALSE.
      REAL, PARAMETER :: excTolerance = 0.2/hartree_to_ev_const
      !-Local Arrays
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      REAL :: ex(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      REAL :: shift(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      REAL :: integrand(gfinp%ne), norm(gfinp%ne)
      COMPLEX, ALLOCATABLE :: imag(:,:,:)
      COMPLEX, ALLOCATABLE :: potmmpmat(:,:,:)


      ALLOCATE(potmmpmat(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, input%jspins))
      ALLOCATE(imag(gfinp%ne,-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const),source=cmplx_0)

      h_loc = 0.0
      DO i_hia = 1, atoms%n_hia

         l     = atoms%lda_u(atoms%n_u+i_hia)%l
         nType = atoms%lda_u(atoms%n_u+i_hia)%atomType

         i_gf = gfinp%hiaElem(i_hia)
         i_elem = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=.TRUE., l_kresolved_int=.TRUE.)
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
            norm = 0.0
            imag = greensfImagPart%applyCutoff(i_elem,i_gf,jspin,.TRUE.)/(3.0-input%jspins)
            DO m = -l, l
               DO mp = -l, l
                  integrand = 0.0
                  DO ie = 1, gfinp%ne
                     integrand(ie) = -1.0/pi_const * ((ie-1) * del+eb) * imag(ie,m,mp)
                     IF(m.EQ.mp) norm(ie) = norm(ie) -1.0/pi_const * imag(ie,m,mp)
                  ENDDO
                  h_loc(m,mp,i_hia,jspin) = trapz(integrand,del,gfinp%ne)
               ENDDO
            ENDDO
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
            potmmpmat = rotMMPmat(v%mmpmat(:,:,i_u,:input%jspins),0.0,-nococonv%beta(nType),-nococonv%alph(nType),l)
         ELSE IF(noco%l_soc) THEN
            potmmpmat = rotMMPmat(v%mmpmat(:,:,i_u,:input%jspins),0.0,-nococonv%theta,-nococonv%phi,l)
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
         IF(input%jspins.EQ.2) THEN
            ex = 0.0
            DO m= -l, l
               DO mp = -l, l
                  ex(m,mp) = h_loc(m,mp,i_hia,1)-h_loc(m,mp,i_hia,2)
               ENDDO
            ENDDO
#ifdef CPP_DEBUG
            WRITE(*,*) "Exchange (eV)"
            WRITE(*,"(7f7.3)") ex(-3:3,-3:3)*hartree_to_ev_const*0.5
#endif
            !------------------------------------------------------------------------------------
            ! If states move close to the cutoff we get some shift in the results on the diagonal
            ! The reason for this is a bit unclear we remove these results and replace them
            ! with either the -m or corresponding opposite spin result (only diagonal)
            !------------------------------------------------------------------------------------
            !Shift m=0 exchange to 0
            IF(l_correctMinus)THEN
               IF(ABS(ex(0,0)).LT.excTolerance) THEN
                  !There is no big discrepancy between m=0 spin up/down => simply eleminate the exchange
                  DO m = -l, l
                     h_loc(m,m,i_hia,1) = h_loc(m,m,i_hia,1) - ex(0,0)/2.0
                     h_loc(m,m,i_hia,2) = h_loc(m,m,i_hia,2) + ex(0,0)/2.0
                  ENDDO
               ELSE
                  !There is a big discrepancy due to numerical problems => Take the spin up part
                  h_loc(0,0,i_hia,2) = h_loc(0,0,i_hia,1)
               ENDIF

               !Recalculate exchange
               ex = 0.0
               DO m= -l, l
                  DO mp = -l, l
                     ex(m,mp) = h_loc(m,mp,i_hia,1)-h_loc(m,mp,i_hia,2)
                  ENDDO
               ENDDO
#ifdef CPP_DEBUG
               WRITE(*,*) "Exchange shifted (eV)"
               WRITE(*,"(7f7.3)") ex(-3:3,-3:3)*hartree_to_ev_const*0.5
#endif
               !Calculate shifts for numerically "troubled" elements
               shift=0.0
               DO m = -l, -1
                  !Normal leftover from SOC+numerical problems
                  IF(ABS(ex(m,m)).LT.excTolerance) CYCLE
                  shift(m,m) = ex(m,m) + ex(-m,-m)
               ENDDO
               h_loc(:,:,i_hia,2) = h_loc(:,:,i_hia,2) + shift

#ifdef CPP_DEBUG
               !Recalculate differences for verification
               ex = 0.0
               DO m= -l, l
                  DO mp = -l, l
                     ex(m,mp) = h_loc(m,mp,i_hia,1)-h_loc(m,mp,i_hia,2)
                  ENDDO
               ENDDO

               WRITE(*,*) "Exchange after Correction (eV)"
               WRITE(*,"(7f7.3)") ex(-3:3,-3:3)*hartree_to_ev_const*0.5
#endif
            ENDIF
        ENDIF

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
   END SUBROUTINE crystal_field

END MODULE m_crystalfield
