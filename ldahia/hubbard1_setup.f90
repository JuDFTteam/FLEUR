MODULE m_hubbard1_setup

   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE hubbard1_setup(iter,atoms,sym,mpi,noco,input,inDen,usdus,pot,gdft,l_runinfleur,l_runhia)

      USE m_types
      USE m_hubbard1_io
      USE m_uj2f
      USE m_constants
      USE m_gfcalc
      USE m_umtx
      USE m_vmmp

      INTEGER,          INTENT(IN)     :: iter !number of iteration 
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_mpi),      INTENT(IN)     :: mpi
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_usdus),    INTENT(IN)     :: usdus
      TYPE(t_potden),   INTENT(IN)     :: inDen
      TYPE(t_potden),   INTENT(INOUT)  :: pot
      TYPE(t_greensf),  INTENT(IN)  :: gdft !green's function in the mt-sphere including the potential correction form dft+hubbard1 
      LOGICAL,          INTENT(IN)     :: l_runinfleur !Determines wether we call the the solver here or run separately
      LOGICAL,          INTENT(IN)     :: l_runhia
      
      INTEGER i_hia,i_gf,n,l,n_occ,ispin,m,matsize,i,iz,N_mu,i_mu,j,io_error
      REAL mu,beta,mu_tmp,a,b,n_tmp,n_max,e_lda_hia
      CHARACTER(len=300) :: cwd
      CHARACTER(len=300) :: path
      CHARACTER(len=300) :: folder
      CHARACTER(len=300) :: message
      TYPE(t_greensf) :: gu,g0

      REAL     f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL     f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      REAL     u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
               -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia)

      COMPLEX  mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX  mmpMat_in(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX  selfen(atoms%n_hia,gdft%nz,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1))
      REAL     e(gdft%nz),dist(input%jspins)
      REAL     n_l(atoms%n_hia,input%jspins)
      LOGICAL  l_exist

      IF(g0%mode.EQ.2) CALL juDFT_error("This energy contour is not supported at the moment for DFT+Hubbard1",calledby="hubbard1_setup")

      
      !the occupation matrix was calculated --> simply calculate potential correction
         !Get the slater integrals from the U and J parameters
      CALL uj2f(input%jspins,atoms%lda_hia(:),atoms%n_hia,f0,f2,f4,f6)
      f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
      f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
      f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
      f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
      CALL umtx(atoms%lda_hia(:),atoms%n_hia,f0(:,1),f2(:,1),f4(:,1),f6(:,1),u)
      IF(ANY(inDen%mmpMat(atoms%n_u+1:atoms%n_hia+atoms%n_u,:,:,:).NE.0.0)) THEN
         mmpMat = inDen%mmpMat(:,:,atoms%n_u+1:atoms%n_hia,:)
         CALL v_mmp(sym,atoms,atoms%lda_hia,atoms%n_hia,input%jspins,.true.,mmpMat,&
         u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),e_lda_hia)
      ELSE !occupation matrix is zero
         !There is nothing to be done yet just set the potential correction to 0
         pot%mmpMat(:,:,atoms%n_u+1:atoms%n_hia+atoms%n_u,:) = CMPLX(0.0,0.0)
      ENDIF

      
      IF(l_runhia.AND.ANY(gdft%gmmpMat(:,:,:,:,:,:,:).NE.0.0)) THEN 
         !The onsite green's function was calculated but the solver 
         !was not yet run
         !--> write out the configuration for the hubbard 1 solver 

         CALL gu%init(input,lmaxU_const,atoms,.true.,noco,nz_in=gdft%nz, e_in=gdft%e,de_in=gdft%de,matsub_in=gdft%nmatsub)
         CALL g0%init(input,lmaxU_const,atoms,.true.,noco,nz_in=gdft%nz, e_in=gdft%e,de_in=gdft%de,matsub_in=gdft%nmatsub)

         !Get the working directory
         CALL get_environment_variable('PWD',cwd)
         !Create a folder where to store the files from the solver
         !Is this applicable on all systems where fleur can be run?
         WRITE(folder,"(A5,I3.3)") "hub1_" , iter
         path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
         CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)))
         !Eliminate the potential V_U from the onsite green's function
         CALL elim_interaction(gdft,g0,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),atoms,usdus,input%jspins)

         DO i_hia = 1, atoms%n_hia
            n = atoms%lda_hia(i_hia)%atomType
            l = atoms%lda_hia(i_hia)%l
            CALL indexgf(atoms,l,n,i_gf)
            matsize = 2*(2*l+1)
            !Create a subdirectory for the atomType and shell
            WRITE(folder,"(A4,I2.2,A2,I1.1)") "atom",n,"_l",l
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
         
            !calculate the occupation of the correlated shell (!!!FROM THE LDA+U GREEN'S FUNCTION)
            CALL occmtx(gdft,i_gf,atoms,sym,input%jspins,input%onsite_beta,mmpMat(:,:,i_hia,:))
            n_l(i_hia,:) = 0.0
            DO ispin = 1, input%jspins
               DO m = -l, l
                  n_l(i_hia,ispin) = n_l(i_hia,ispin) + mmpMat(m,m,i_hia,ispin)
               ENDDO
            ENDDO
            WRITE(*,*) "OCCUPATION: ", SUM(n_l(i_hia,:))
            !IF(ALL(pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:).EQ.0.0)) THEN
            !   n_l(i_hia,1) = 7.0
            !   n_l(i_hia,2) = 0.0
            !ENDIF 

            !calculate the chemical potential
            CALL mudc(atoms%lda_hia(i_hia)%U,atoms%lda_hia(i_hia)%J,l,n_l(i_hia,:),atoms%lda_hia(i_hia)%l_amf,mu,input%jspins)

            !minimum and maximum number of electrons considered in the solver (+- 2 from nearest integer)
            
            n_occ = ANINT(SUM(n_l(i_hia,:)))
            n_l(i_hia,1) = SUM(n_l(i_hia,:))
            !CHECK wether the hubbard 1 solver was run:
            INQUIRE(file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/" // "se.atom",exist=l_exist)

            IF(l_exist) THEN
               CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,1:g0%nz-g0%nmatsub,1:2*(2*l+1),1:2*(2*l+1)),g0%nz-g0%nmatsub,2*(2*l+1),e(:),.false.)
               CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,g0%nz-g0%nmatsub:g0%nz,1:2*(2*l+1),1:2*(2*l+1)),g0%nmatsub,2*(2*l+1),e(:),.true.)
            ELSE
               CALL write_hubbard1_input(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),&
                                          0.20,-0.0001,MAX(1,n_occ-2),MIN(2*(2*l+1),n_occ+2),input%onsite_beta,mu,&
                                          gdft%nz-gdft%nmatsub,gdft%nmatsub,REAL(g0%e(1)),REAL(g0%e(g0%nz-g0%nmatsub)),AIMAG(g0%e(1)))

               OPEN(unit=1337,file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/" // "contour.dat",status="replace",action="write")

               DO iz = 1, gdft%nz
                  WRITE(1337,"(2f14.8)") REAL(gdft%e(iz))*hartree_to_ev_const, AIMAG(gdft%e(iz))*hartree_to_ev_const
               ENDDO

               CLOSE(unit= 1337)
               
               !There has to be a better solution
               !Maybe use CALL System() to start the solver from here
               !EXPERIMENTAL:
               IF(l_runinfleur) THEN
                  CALL timestart("Hubbard 1 solver")
                  CALL CHDIR(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/")
                  WRITE(*,"(A)") "Running Hubbard 1 solver"
                  CALL SYSTEM("/home/henning/GIT/hub2new4sp/eigen > out")
                  CALL SYSTEM("/home/henning/GIT/hub2new4sp/selfen > out") 
                  WRITE(*,"(A)") "Hubbard 1 solver finished" 
                  CALL timestop("Hubbard 1 solver")
                  CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,1:g0%nz-g0%nmatsub,1:2*(2*l+1),1:2*(2*l+1)),g0%nz-g0%nmatsub,2*(2*l+1),e(:),.false.)
                  IF(g0%nmatsub.GT.0) CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,g0%nz-g0%nmatsub:g0%nz,1:2*(2*l+1),1:2*(2*l+1)),g0%nmatsub,2*(2*l+1),e(:),.true.)
               ENDIF
            ENDIF
         ENDDO
         CALL CHDIR(TRIM(ADJUSTL(cwd)))

         IF(l_exist.OR.l_runinfleur) THEN
            CALL add_selfen(g0,gu,selfen,atoms,sym,input%onsite_beta,input%jspins,n_l(:,1),mmpMat)
            ! calculate potential matrix and total energy correction
            CALL v_mmp(sym,atoms,atoms%lda_hia,atoms%n_hia,input%jspins,.true.,mmpMat,&
                  u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),e_lda_hia)
            !
            ! Output the density of states from the two green's functions
            !
            DO i_hia = 1, atoms%n_hia
               n = atoms%lda_hia(i_hia)%atomType
               l = atoms%lda_hia(i_hia)%l
               WRITE(folder,"(A5,I3.3)") "hub1_" , iter
               path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
               WRITE(folder,"(A4,I2.2,A2,I1.1)") "atom",n,"_l",l
               CALL CHDIR(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/")
               CALL indexgf(atoms,l,n,i_gf)
               CALL ldosmtx("gdft",gdft,i_gf,atoms,sym,input%jspins)
               CALL ldosmtx("g0",g0,i_gf,atoms,sym,input%jspins)
               CALL ldosmtx("g",gu,i_gf,atoms,sym,input%jspins)
            ENDDO
            CALL CHDIR(TRIM(ADJUSTL(cwd)))
            INQUIRE(file="n_mmpmat_hubbard1",exist = l_exist)
            IF(l_exist) THEN
               OPEN(unit=1337,file="n_mmpmat_hubbard1",status="old",action="read",iostat=io_error)
               IF(io_error.NE.0) CALL juDFT_error("IO-error in density matrix",calledby="hubbard1_setup")
               READ(1337,"(7f14.8)") mmpMat_in(:,:,:,:)
               CLOSE(unit=1337)
               CALL n_mmp_dist(mmpMat_in,mmpMat,atoms,dist,input%jspins)
               inquire(file="n_mmp_dist", exist=l_exist)
               if (l_exist) then
                 open(12, file="n_mmp_dist", status="old", position="append", action="write")
               else
                 open(12, file="n_mmp_dist", status="new", action="write")
               end if
               write(12, "(f14.8)") NORM2(dist)
               close(12)
            ENDIF 
  
            !Write out the density matrix (because i don't wanna change inDen to intent(inout)
            OPEN(unit=1337,file="n_mmpmat_hubbard1",status="replace",action="write",iostat=io_error)
            IF(io_error.NE.0) CALL juDFT_error("IO-error in density matrix",calledby="hubbard1_setup")
            WRITE(1337,"(7f14.8)") mmpMat(:,:,:,:)
            CLOSE(unit=1337)
            IF(l_exist) THEN
               IF(NORM2(dist).LT.0.001) THEN
                  CALL juDFT_END("Density matrix has converged")
               ENDIF
            ENDIF
         ELSE
            WRITE(message,9010) iter
            CALL juDFT_END(message)
9010        FORMAT("Hubbard 1 input written into hub1_",I3.3)
         ENDIF
      END IF

      WRITE(*,"(7f14.8)") REAL(pot%mmpMat(:,:,:,:))
     
   END SUBROUTINE hubbard1_setup

   SUBROUTINE mudc(U,J,l,n,l_amf,mu,jspins)

      REAL,    INTENT(IN)  :: U
      REAL,    INTENT(IN)  :: J 
      INTEGER, INTENT(IN)  :: l
      REAL,    INTENT(IN)  :: n(jspins)
      LOGICAL, INTENT(IN)  :: l_amf
      REAL,    INTENT(OUT) :: mu
      INTEGER, INTENT(IN)  :: jspins

      REAL vdcfll1,vdcfll2
      REAL vdcamf1,vdcamf2
      REAL nup,ndn

      IF(jspins.EQ.2) THEN
         nup = n(1)
         ndn = n(2)
      ELSE
         nup = 0.5 * n(1)
         ndn = nup
      ENDIF

      vdcfll1= u*(nup+ndn -0.5) - j*(nup-0.5)
      vdcfll2= u*(nup+ndn -0.5) - j*(ndn-0.5)
      vdcamf1= u*ndn+2*l/(2*l+1)*(u-j)*nup
      vdcamf2= u*nup+2*l/(2*l+1)*(u-j)*ndn
      WRITE(*,"(A)") 'Double counting chemical potential:'
      WRITE(*,9010) 'FLL: ','spin-up','spin-dn','(up+dn)/2','up-dn'
      WRITE(*,9020) vdcfll1,vdcfll2,(vdcfll1+vdcfll2)/2,vdcfll1-vdcfll2
      WRITE(*,9010) 'AMF: ','spin-up','spin-dn','(up+dn)/2','up-dn'
      WRITE(*,9020)  vdcamf1,vdcamf2,(vdcamf1+vdcamf2)/2,vdcamf1-vdcamf2

      IF(l_amf) THEN
         WRITE(*,"(A)") "Using the around-mean-field limit"
         mu = (vdcamf1+vdcamf2)/2
      ELSE
      WRITE(*,"(A)") "Using the fully-localized limit"
         mu = (vdcfll1+vdcfll2)/2
      ENDIF 
      WRITE(*,9030) mu

9010  FORMAT(TR3,A4,TR1,A7,TR3,A7,TR3,A9,TR3,A5)
9020  FORMAT(TR7,f8.4,TR2,f8.4,TR2,f8.4,TR4,f8.4)
9030  FORMAT(TR3,"mu = ",f7.4)
   END SUBROUTINE

   SUBROUTINE add_selfen(g,gp,selfen,atoms,sym,beta,jspins,n_occ,mmpMat)

      !Calculates the interacting Green's function for the mt-sphere with
      !
      ! (G)^-1 = (G_0)^-1 - mu 1 - selfen
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

      TYPE(t_greensf),  INTENT(IN)     :: g
      TYPE(t_greensf),  INTENT(INOUT)  :: gp
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      COMPLEX,          INTENT(IN)     :: selfen(atoms%n_hia,g%nz,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1))
      INTEGER,          INTENT(IN)     :: jspins
      REAL,             INTENT(IN)     :: beta
      REAL,             INTENT(IN)     :: n_occ(atoms%n_hia)
      COMPLEX,          INTENT(OUT)    :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,jspins)

      INTEGER i_hia,l,nType,i_gf,ns,ispin,m,iz,ipm
      CHARACTER(len=6) app

      REAL mu_a,mu_b,mu_step,mu_max,n_max
      REAL mu,n

      TYPE(t_mat) :: gmat,vmat

      !Interval where we expect the correct mu
      mu_a = -2.0
      mu_b = 2.0
      mu_step = 0.05
      mu_max = 0.0
      n_max = 0.0

      DO i_hia = 1, atoms%n_hia
         l = atoms%lda_hia(i_hia)%l
         nType = atoms%lda_hia(i_hia)%atomType
         ns = 2*l+1
         CALL indexgf(atoms,l,nType,i_gf)
         !intialize the matrices
         CALL gmat%init(.false.,2*ns,2*ns)
         CALL vmat%init(.false.,2*ns,2*ns)
         !Search for the maximum of occupation
         OPEN(unit=1337,file="mu",status="replace",action="write")
         mu = mu_a
         DO WHILE(mu.LE.mu_b)
            mu = mu + mu_step
            DO iz = 1, g%nz
               vmat%data_c = selfen(i_hia,iz,1:2*ns,1:2*ns)
               DO ipm = 1, 2
                  CALL to_tmat(gmat,g%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,l)
                  CALL add_pot(gmat,vmat,mu,(ipm.EQ.1))
                  CALL to_g(gmat,gp%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,l)
               ENDDO
            ENDDO
            CALL occmtx(gp,i_gf,atoms,sym,jspins,beta,mmpMat(:,:,i_hia,:))
            !Calculate the trace
            n = 0.0
            DO ispin = 1, jspins
               DO m = -l, l
                  n = n + mmpMat(m,m,i_hia,ispin)
               ENDDO
            ENDDO
            WRITE(1337,*) mu,n
            IF(n.GT.n_max) THEN
               mu_max = mu
               n_max  = n
            ENDIF
         ENDDO
         CLOSE(1337)
         !Set up the interval for the bisection method (mu_max,mu_b)
         mu_a = mu_max
         DO 
            mu = (mu_a + mu_b)/2.0
            DO iz = 1, g%nz
               vmat%data_c = selfen(i_hia,iz,1:2*ns,1:2*ns)
               DO ipm = 1, 2
                  CALL to_tmat(gmat,g%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,l)
                  CALL add_pot(gmat,vmat,mu,(ipm.EQ.1))
                  CALL to_g(gmat,gp%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,l)
               ENDDO
            ENDDO
            CALL occmtx(gp,i_gf,atoms,sym,jspins,beta,mmpMat(:,:,i_hia,:))
            !Calculate the trace
            n = 0.0
            DO ispin = 1, jspins
               DO m = -l, l
                  n = n + mmpMat(m,m,i_hia,ispin)
               ENDDO
            ENDDO
            IF(ABS(n-n_occ(i_hia)).LT.0.001.OR.ABS((mu_b - mu_a)/2.0).LT.0.00001) THEN
               !We found the chemical potential to within the desired accuracy
               !TODO: Write to output file
               WRITE(*,*) "Calculated mu: ", mu
               WRITE(*,*) "OCCUPATION: ", n
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


   SUBROUTINE elim_interaction(g,gp,v,atoms,usdus,jspins)

      !Calculates the new non-interacting green's function for 
      !DFT+Hubbard 1 by eliminating the potential correction with
      !
      ! (G_0)^-1 = (G_U)^-1 + V_U
      !
      !where G_U is the green's function calculated from the KS-eigenstates with V_U

      USE m_types
      USE m_gfcalc
      USE m_constants

      TYPE(t_greensf),     INTENT(IN)     :: g
      TYPE(t_greensf),     INTENT(INOUT)  :: gp
      TYPE(t_atoms),       INTENT(IN)     :: atoms 
      COMPLEX,             INTENT(IN)     :: v(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,jspins)
      TYPE(t_usdus),       INTENT(IN)     :: usdus
      INTEGER,             INTENT(IN)     :: jspins

      INTEGER i_hia,iz,ipm,i_gf,l,n,ns

      TYPE(t_mat) :: gmat,vmat

      DO i_hia = 1, atoms%n_hia
         l = atoms%lda_hia(i_hia)%l
         n = atoms%lda_hia(i_hia)%atomType
         ns = 2*l+1
         CALL indexgf(atoms,l,n,i_gf)
         !intialize the matrices
         CALL gmat%init(.false.,2*ns,2*ns)
         CALL vmat%init(.false.,2*ns,2*ns)
         CALL to_tmat(vmat,-v(:,:,i_hia,:),jspins,l)
         WRITE(*,"(14f10.5)") REAL(vmat%data_c)
         DO iz = 1, g%nz
            DO ipm = 1, 2
               CALL to_tmat(gmat,g%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,l)
               CALL add_pot(gmat,vmat,0.0,(ipm.EQ.1))
               CALL to_g(gmat,gp%gmmpMat(1,iz,i_gf,:,:,:,ipm),jspins,l)
            ENDDO
         ENDDO
         CALL gmat%free()
         CALL vmat%free()
      ENDDO

   END SUBROUTINE elim_interaction

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

   SUBROUTINE n_mmp_dist(n_mmp_in,n_mmp_out,atoms,dist,jspins)

      USE m_types
      USE m_constants

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      COMPLEX,             INTENT(IN)  :: n_mmp_in(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,jspins)
      COMPLEX,             INTENT(IN)  :: n_mmp_out(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,jspins)
      REAL,                INTENT(OUT) :: dist(jspins)
      INTEGER,             INTENT(IN)  :: jspins

      INTEGER ispin,i_hia,j,k
      
      !Calculates the distance for two density matrices
      DO ispin = 1, jspins
         dist(ispin) = 0.0
         DO i_hia = 1, atoms%n_hia
            DO j = -3,3
               DO k = -3,3
                  dist(ispin) = dist(ispin) + ABS(n_mmp_out(k,j,i_hia,1) - n_mmp_in(k,j,i_hia,1))
               END DO
            END DO
         END DO
         WRITE(*,9010) ispin, dist(ispin)
9010     FORMAT("n_mmp_distance for spin ",I1.1,": ",f12.6)
      ENDDO
   
   END SUBROUTINE n_mmp_dist

END MODULE m_hubbard1_setup