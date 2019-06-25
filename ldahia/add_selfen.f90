MODULE m_add_selfen

   LOGICAL, PARAMETER :: l_debug = .TRUE.

   CONTAINS

   SUBROUTINE add_selfen(g,gp,selfen,atoms,hub1,sym,input,ef,n_occ,mu_dc,vmmp,mmpMat)

      !Calculates the interacting Green's function for the mt-sphere with
      !
      ! (G)^-1 = (G_0)^-1 - mu 1 + V_dc - selfen
      !
      !The term mu * unity is there to ensure that the number of particles 
      !doesn't change and is determined by a two-step process
      !The occupation as a function of mu has a peak in the region where 
      !something is inside the energy interval between e_bot adn e_fermi
      !To determine where we have the same number of particles we first 
      !search for the maximum occupation
      !Then the desired chemical potential is found with the bisection method 
      !to the right of the maximum
      
      USE m_types
      USE m_constants
      USE m_gfcalc

      IMPLICIT NONE

      TYPE(t_greensf),  INTENT(IN)     :: g
      TYPE(t_greensf),  INTENT(INOUT)  :: gp
      TYPE(t_hub1ham),  INTENT(IN)     :: hub1
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_input),    INTENT(IN)     :: input
      COMPLEX,          INTENT(IN)     :: selfen(hub1%n_hia,g%nz,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1))
      REAL,             INTENT(IN)     :: ef
      REAL,             INTENT(IN)     :: n_occ(hub1%n_hia)
      REAL,             INTENT(IN)     :: mu_dc
      COMPLEX,          INTENT(IN)     :: vmmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,input%jspins)
      COMPLEX,          INTENT(OUT)    :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,input%jspins)

      INTEGER i_hia,l,nType,i_gf,ns,ispin,m,iz,ipm
      CHARACTER(len=6) app

      REAL mu_a,mu_b,mu_step,mu_max,n_max
      REAL mu,n
      LOGICAL l_mperp

      TYPE(t_mat) :: gmat,vmat

      !replace with noco%l_mperp
      l_mperp = .true.
      !Interval where we expect the correct mu
      mu_a = -4.0
      mu_b = 4.0
      mu_step = 0.05
      mu_max = 0.0
      n_max = 0.0

      DO i_hia = 1, hub1%n_hia
         l = hub1%lda_u(i_hia)%l
         nType = hub1%lda_u(i_hia)%atomType
         ns = 2*l+1
         CALL indexgf(atoms,l,nType,i_gf)
         !intialize the matrices
         CALL gmat%init(.false.,2*ns,2*ns)
         CALL vmat%init(.false.,2*ns,2*ns)
         !Search for the maximum of occupation
         IF(l_debug) OPEN(unit=1337,file="mu",status="replace",action="write")
         mu = mu_a
         DO WHILE(mu.LE.mu_b)
            mu = mu + mu_step
            DO iz = 1, g%nz
               vmat%data_c = selfen(i_hia,iz,1:2*ns,1:2*ns)
               IF(.NOT.l_mperp) THEN
                  !Dismiss spin-off-diagonal elements
                  vmat%data_c(1:ns,ns+1:2*ns) = 0.0
                  vmat%data_c(ns+1:2*ns,1:ns) = 0.0
               ENDIF
               DO ipm = 1, 2
                  CALL to_tmat(gmat,g%gmmpMat(iz,i_gf,:,:,:,ipm),input%jspins,2,l)
                  CALL add_pot(gmat,vmat,mu-mu_dc,(ipm.EQ.1))
                  CALL to_g(gmat,gp%gmmpMat(iz,i_gf,:,:,:,ipm),2,input%jspins,l)
               ENDDO
            ENDDO
            CALL occmtx(gp,i_gf,atoms,sym,input,ef,mmpMat(:,:,i_hia,:))
            !Calculate the trace
            n = 0.0
            DO ispin = 1, input%jspins
               DO m = -l, l
                  n = n + mmpMat(m,m,i_hia,ispin)
               ENDDO
            ENDDO
            IF(l_debug) WRITE(1337,"(2f15.8)") mu,n
            IF(n.GT.n_max) THEN
               mu_max = mu
               n_max  = n
            ENDIF
         ENDDO
         IF(l_debug) CLOSE(1337)

         !Sanity check for the maximum occupation
         IF(n_max.GT.2*ns+0.5) THEN
            !These oscillations seem to emerge when the lorentzian smoothing is done inadequately
            CALL juDFT_error("Something went wrong with the addition of the selfenergy",calledby="add_selfen")
         ENDIF

         !Set up the interval for the bisection method (mu_max,mu_b)
         mu_a = mu_max
         DO 
            mu = (mu_a + mu_b)/2.0
            DO iz = 1, g%nz
               vmat%data_c = selfen(i_hia,iz,1:2*ns,1:2*ns)
               IF(.NOT.l_mperp) THEN
                  !Dismiss spin-off-diagonal elements
                  vmat%data_c(1:ns,ns+1:2*ns) = 0.0
                  vmat%data_c(ns+1:2*ns,1:ns) = 0.0
               ENDIF
               DO ipm = 1, 2
                  CALL to_tmat(gmat,g%gmmpMat(iz,i_gf,:,:,:,ipm),input%jspins,2,l)
                  CALL add_pot(gmat,vmat,mu-mu_dc,(ipm.EQ.1))
                  CALL to_g(gmat,gp%gmmpMat(iz,i_gf,:,:,:,ipm),2,input%jspins,l)
               ENDDO
            ENDDO
            CALL occmtx(gp,i_gf,atoms,sym,input,ef,mmpMat(:,:,i_hia,:))
            !Calculate the trace
            n = 0.0
            DO ispin = 1, input%jspins
               DO m = -l, l
                  n = n + mmpMat(m,m,i_hia,ispin)
               ENDDO
            ENDDO
            IF(ABS(n-n_occ(i_hia)).LT.0.001.OR.ABS((mu_b - mu_a)/2.0).LT.0.00001) THEN
               !We found the chemical potential to within the desired accuracy
               WRITE(6,"(A)") "Calculated mu to match Self-energy to DFT-GF"
               WRITE(6,"(TR3,A4,f8.4)") "mu = ", mu
               EXIT
            ELSE IF((n - n_occ(i_hia)).GT.0) THEN
               !The occupation is to small --> choose the left interval
               mu_a = mu
            ELSE IF((n - n_occ(i_hia)).LT.0) THEN
               !The occupation is to big --> choose the right interval
               mu_b = mu
            ENDIF
         ENDDO
         CALL gmat%free()
         CALL vmat%free()
      ENDDO

   END SUBROUTINE add_selfen

   SUBROUTINE add_pot(gmat,vmat,mu,l_upper)

      USE m_types

      TYPE(t_mat),      INTENT(INOUT)  :: gmat
      TYPE(t_mat),      INTENT(IN)     :: vmat
      REAL,             INTENT(IN)     :: mu
      LOGICAL,          INTENT(IN)     :: l_upper !Are we in the upper half of the complex plane

      INTEGER i,j
      
      CALL gmat%inverse()
      DO i = 1, gmat%matsize1
         DO j = 1, gmat%matsize1
            IF(l_upper) THEN
               gmat%data_c(i,j) = gmat%data_c(i,j) - vmat%data_c(i,j)
            ELSE
               gmat%data_c(i,j) = gmat%data_c(i,j) - conjg(vmat%data_c(i,j))
            ENDIF
            IF(i.EQ.j) gmat%data_c(i,i) = gmat%data_c(i,i) - mu
         ENDDO
      ENDDO
      CALL gmat%inverse()
    
   END SUBROUTINE add_pot

END MODULE m_add_selfen