MODULE m_add_selfen

   USE m_types
   USE m_types_selfen
   USE m_constants
   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE add_selfen(g0,selfen,gfinp,input,atoms,noco,nococonv,occDFT,g,mmpMat)

      !Calculates the interacting Green's function for the mt-sphere with
      !
      ! (G)^-1 = (G_0)^-1 - mu 1 + V_dc - selfen + V_U
      !
      !The term mu * unity is there to ensure that the number of particles
      !doesn't change and is determined by a two-step process
      !The occupation as a function of mu has a peak in the region where
      !something is inside the energy interval between e_bot and e_fermi
      !To determine where we have the same number of particles we first
      !search for the maximum occupation
      !Then the desired chemical potential is found with the bisection method
      !to the right of the maximum
      !TODO: Parallelization (OMP over chemical potentials in first loop??)

      TYPE(t_greensf),  INTENT(IN)     :: g0
      TYPE(t_gfinp),    INTENT(IN)     :: gfinp
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_nococonv), INTENT(IN)     :: nococonv
      REAL,             INTENT(IN)     :: occDFT(:)
      TYPE(t_selfen),   INTENT(INOUT)  :: selfen
      TYPE(t_greensf),  INTENT(INOUT)  :: g
      COMPLEX,          INTENT(INOUT)  :: mmpMat(-lmaxU_const:,-lmaxU_const:,:)

      INTEGER :: l,ispin,m,mp,iMatch
      REAL    :: mu_a,mu_b,mu_step
      REAL    :: mu,nocc,nTarget,muMax,nMax
      LOGICAL :: l_fullMatch,l_invalidElements

      !Are we matching the spin polarized self-energy with one chemical potential
      l_fullMatch = SIZE(selfen%muMatch).EQ.1 .AND. input%jspins.EQ.2

      l = g0%elem%l

      !Search for the maximum of occupation
      DO iMatch = 1, SIZE(selfen%muMatch)

         !Target occupation
         nTarget = MERGE(SUM(occDFT(:)),occDFT(iMatch),l_fullMatch)

         !Interval where we expect the correct mu
         mu_a = -2.0
         mu_b = 1.5
         mu_step = 0.01

         !----------------------------------------------------
         ! Scan the given Interval for the maximum Occupation
         !----------------------------------------------------
         mu = mu_a
         muMax = 0.0
         nMax  = 0.0
         DO WHILE(mu.LE.mu_b)

            mu = mu + mu_step

            CALL getOccupationMtx(g0,gfinp,input,atoms,noco,nococonv,selfen,mu,l_fullMatch,iMatch,&
                                  g,mmpMat,nocc,l_invalidElements)

            IF(nocc.GT.nMax) THEN
               muMax = mu
               nMax  = nocc
            ENDIF

#ifdef CPP_DEBUG
            WRITE(1337,'(2f15.8)') mu,nocc
#endif

         ENDDO

         !Sanity check for the maximum occupation
         IF(nMax-2*(2*l+1).GT.1.0) THEN
            !These oscillations seem to emerge when the lorentzian smoothing is done inadequately
            CALL juDFT_error("Something went wrong with the addition of the selfenergy: nMax>>2*(2*l+1)",&
                              calledby="add_selfen")
         ELSE IF(nMax-nTarget.LT.-0.1) THEN
            CALL juDFT_error("Something went wrong with the addition of the selfenergy: nMax<nTarget",&
                              calledby="add_selfen")
         ENDIF

         !------------------------------------------------------------------
         ! Find the matching chemical potential on the right of the maximum
         !------------------------------------------------------------------
         !Set up the interval for the bisection method (mu_max,mu_b)
         mu_a = muMax
         DO WHILE (ABS(nocc-nTarget).GT.1e-8.AND.ABS((mu_b - mu_a)/2.0).GT.1e-8)
            mu = (mu_a + mu_b)/2.0

            CALL getOccupationMtx(g0,gfinp,input,atoms,noco,nococonv,selfen,mu,l_fullMatch,iMatch,&
                                  g,mmpMat,nocc,l_invalidElements)

            IF((nocc - nTarget).GT.0.0) THEN
               !The occupation is to big --> choose the right interval
               mu_a = mu
            ELSE IF((nocc - nTarget).LT.0.0) THEN
               !The occupation is to small --> choose the left interval
               mu_b = mu
            ENDIF
         ENDDO
         selfen%muMatch(iMatch) = mu
         !----------------------------------------------------
         ! Check if the final mmpMat contains invalid elements
         !----------------------------------------------------
         IF(l_invalidElements) THEN
            CALL juDFT_error("Invalid Element/occupation in final density matrix",calledby="add_selfen")
         ENDIF
      ENDDO

      !Test throw out elements smaller than 1e-4
      DO ispin = 1, MERGE(3,input%jspins,gfinp%l_mperp)
         DO m = -l, l
            DO mp =-l, l
               IF(ABS(mmpMat(m,mp,ispin)).LT.1e-4) mmpMat(m,mp,ispin) = cmplx_0
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE add_selfen

   SUBROUTINE getOccupationMtx(g0,gfinp,input,atoms,noco,nococonv,selfen,mu,l_fullMatch,iMatch,&
                               g,mmpMat,nocc,l_invalidElements)


      TYPE(t_greensf),     INTENT(IN)    :: g0           !DFT Green's Function
      TYPE(t_gfinp),       INTENT(IN)    :: gfinp
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_noco),        INTENT(IN)    :: noco
      TYPE(t_nococonv),    INTENT(IN)    :: nococonv
      TYPE(t_selfen),      INTENT(IN)    :: selfen       !Atomic self-energy (with removed LDA+U potential)
      REAL,                INTENT(IN)    :: mu           !chemical potential shift
      LOGICAL,             INTENT(IN)    :: l_fullMatch  !Are spins matched individually?
      INTEGER,             INTENT(IN)    :: iMatch       !Index for current matching
      TYPE(t_greensf),     INTENT(INOUT) :: g            !Impurity Green's Function
      COMPLEX,             INTENT(INOUT) :: mmpMat(-lmaxU_const:,-lmaxU_const:,:) !Occupation matrix of g
      REAL,                INTENT(INOUT) :: nocc         !trace of the occupation matrix
      LOGICAL,             INTENT(INOUT) :: l_invalidElements !Are there invalid elements in the resulting Occupation matrix

      INTEGER :: l,ns,matsize,start,end
      INTEGER :: ipm,iz,ispin,m

      TYPE(t_mat) :: vmat,gmat

      l = g0%elem%l
      ns = 2*l+1
      matsize = ns*MERGE(2,1,l_fullMatch)
      CALL vmat%init(.false.,matsize,matsize)
      CALL gmat%init(.false.,matsize,matsize)
      IF(iMatch>1) CALL g%reset()

      !Select the correct section from the selfenergy
      start = MERGE(1,1+(iMatch-1)*ns,l_fullMatch)
      end   = MERGE(2*ns,iMatch*ns,l_fullMatch)
      DO ipm = 1, 2
         DO iz = 1, g0%contour%nz
            !Read selfenergy
            vmat%data_c = selfen%data(start:end,start:end,iz,ipm)
            IF(.NOT.gfinp%l_mperp.AND.l_fullMatch) THEN
               !Dismiss spin-off-diagonal elements
               vmat%data_c(1:ns,ns+1:2*ns) = cmplx_0
               vmat%data_c(ns+1:2*ns,1:ns) = cmplx_0
            ENDIF

            !Read in the DFT-Green's Function at the energy point
            IF(l_fullMatch) THEN
               CALL g0%getFullMatrix(atoms,iz,ipm.EQ.2,gmat)
            ELSE
               CALL g0%get(atoms,iz,ipm.EQ.2,iMatch,gmat)
            ENDIF

            !----------------------------------------------------
            !Solve the Dyson equation at the current energy point
            !----------------------------------------------------
            CALL add_pot(gmat,vmat,mu)

            !Set the Impurity-Green's Function at the energy point
            IF(l_fullMatch) THEN
               CALL g%set(iz,ipm.EQ.2,gmat)
            ELSE
               CALL g%set(iz,ipm.EQ.2,gmat,spin=iMatch)
            ENDIF

         ENDDO
      ENDDO

      !Get the occupation matrix
      IF(l_fullMatch) THEN
         mmpMat = g%occupationMatrix(gfinp,input,atoms,noco,nococonv,check=.TRUE.,occError=l_invalidElements)
      ELSE
         mmpMat(:,:,iMatch) = g%occupationMatrix(iMatch,gfinp,input,atoms,noco,nococonv,check=.TRUE.,occError=l_invalidElements)
      ENDIF

      !Compute the trace
      nocc = 0.0
      DO ispin = MERGE(1,iMatch,l_fullMatch), MERGE(input%jspins,iMatch,l_fullMatch)
         DO m = -l, l
            nocc = nocc + REAL(mmpMat(m,m,ispin))
         ENDDO
      ENDDO


   END SUBROUTINE getOccupationMtx

   SUBROUTINE add_pot(gmat,vmat,mu)

      TYPE(t_mat),      INTENT(INOUT)  :: gmat
      TYPE(t_mat),      INTENT(IN)     :: vmat
      REAL,             INTENT(IN)     :: mu

      INTEGER i

      IF(vmat%matsize1.NE.gmat%matsize1.OR.vmat%matsize2.NE.gmat%matsize2) &
         CALL juDFT_error("vmat & gmat dimension do not match",hint="This is a bug in FLEUR, please report",calledby="add_pot")

      CALL gmat%inverse()
      gmat%data_c = gmat%data_c - vmat%data_c
      DO i = 1, gmat%matsize1
         gmat%data_c(i,i) = gmat%data_c(i,i) - mu
      ENDDO
      CALL gmat%inverse()

   END SUBROUTINE add_pot

END MODULE m_add_selfen
