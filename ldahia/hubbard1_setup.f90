MODULE m_hubbard1_setup

   USE m_juDFT
   USE m_types
   USE m_hubbard1_io
   USE m_uj2f
   USE m_constants
   USE m_gfcalc

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE hubbard1_setup(iter,atoms,sym,mpi,input,inDen,pot,gOnsite,l_runinfleur)

      INTEGER,          INTENT(IN)     :: iter !number of iteration 
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_mpi),      INTENT(IN)     :: mpi
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_potden),   INTENT(IN)     :: inDen
      TYPE(t_potden),   INTENT(INOUT)  :: pot
      TYPE(t_greensf),  INTENT(IN)     :: gOnsite
      LOGICAL,          INTENT(IN)     :: l_runinfleur !Determines wether we call the the solver here or run separately
      
      INTEGER i_hia,i_gf,n,l,n_occ,ispin,m
      REAL mu,beta
      CHARACTER(len=300) :: cwd
      CHARACTER(len=300) :: path
      CHARACTER(len=300) :: folder

      REAL     f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL     f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      COMPLEX  mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins)
      REAL     n_l(input%jspins)

      beta = 10.0
      !Check wether the hubbard 1 solver was run for this iteration
      !Then there should be a occupation matrix file tagged with the current iteration
      IF(ANY(inDen%mmpMat(atoms%n_u+1:atoms%n_hia+atoms%n_u,:,:,:).NE.0.0)) THEN
         !the occupation matrix was calculated
      ELSE IF(ANY(gOnsite%gmmpMat(:,:,:,:,:,:,:).NE.0.0)) THEN 
         !The onsite green's function was calculated but the solver 
         !was not yet run
         !--> write out the configuration for the hubbard 1 solver 

         !Get the working directory
         CALL get_environment_variable('PWD',cwd)
         !Create a folder where to store the files from the solver
         !Is this applicable on all systems where fleur can be run?
         WRITE(folder,"(A5,I3.3)") "hub1_" , iter
         path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
         CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)))

         !Get the slater integrals from the U and J parameters
         CALL uj2f(input%jspins,atoms%lda_hia(:),atoms%n_hia,f0,f2,f4,f6)
         f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
         f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
         f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
         f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2

         DO i_hia = 1, atoms%n_hia

            n = atoms%lda_hia(i_hia)%atomType
            l = atoms%lda_hia(i_hia)%l
            !Create a subdirectory for the atomType and shell
            WRITE(folder,"(A4,I2.2,A2,I1.1)") "atom",n,"_l",l
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
            
            !calculate the occupation of the correlated shell
            CALL indexgf(atoms,l,n,i_gf)
            CALL occmtx(gOnsite,i_gf,atoms,sym,input%jspins,mmpMat)
            n_l = 0.0
            DO ispin = 1, input%jspins
               DO m = -l, l
                  n_l(ispin) = n_l(ispin) + mmpMat(m,m,ispin)
               ENDDO
            ENDDO

            !calculate the chemical potential
            CALL mudc(atoms%lda_hia(i_hia)%U,atoms%lda_hia(i_hia)%J,l,n_l(1),n_l(2),atoms%lda_hia(i_hia)%l_amf,mu)

            !minimum and maximum number of electrons considered in the solver (+- 2 from nearest integer)
            n_occ = ANINT(n_l(1)+n_l(2))

            CALL write_hubbard1_input(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),&
                                       0.16,-0.00015,MAX(1,n_occ-2),MIN(2*(2*l+1),n_occ+2),beta,mu)
            !There has to be a better solution
            !Maybe use CALL System() to start the solver from here
            !EXPERIMENTAL:
            IF(l_runinfleur) THEN
               CALL CHDIR(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/")
               CALL SYSTEM("/home/henning/GIT/hub2new4sp/eigen")
               CALL SYSTEM("/home/henning/GIT/hub2new4sp/dos")
               CALL SYSTEM("/home/henning/GIT/hub2new4sp/selfen")
               CALL SYSTEM("/home/henning/GIT/hub2new4sp/angmom")    
            ELSE
               !If we do this we need to write out the green's function in some form
               CALL juDFT_error("The input configurations for the hubbard 1 solver have been written: Please execute the solver now",calledby="hubbard1_setup")
            ENDIF
         ENDDO
         CALL CHDIR(TRIM(ADJUSTL(cwd)))
      ELSE
         IF(mpi%irank.EQ.0) WRITE(*,*) "no onsite green's function found ... skipping LDA+Hubbard 1"
         pot%mmpMat(atoms%n_u+1:atoms%n_hia+atoms%n_u,:,:,:) = CMPLX(0.0,0.0)
      ENDIF

   END SUBROUTINE hubbard1_setup

   SUBROUTINE mudc(U,J,l,nup,ndn,l_amf,mu)

      REAL,    INTENT(IN)  :: U
      REAL,    INTENT(IN)  :: J 
      INTEGER, INTENT(IN)  :: l
      REAL,    INTENT(IN)  :: nup
      REAL,    INTENT(IN)  :: ndn
      LOGICAL, INTENT(IN)  :: l_amf
      REAL,    INTENT(OUT) :: mu

      REAL vdcfll1,vdcfll2
      REAL vdcamf1,vdcamf2

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

   SUBROUTINE indexgf(atoms,l,n,ind)

      !Find the index of the greens function associated with this l,n 

      USE m_juDFT
      USE m_types

      !Finds the corresponding entry in gmmpMat for given atomType and l

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      INTEGER,             INTENT(IN)  :: l,n
      INTEGER,             INTENT(OUT) :: ind

      ind = 0
      DO 
         ind = ind + 1
         IF(atoms%onsiteGF(ind)%atomType.EQ.n.AND.atoms%onsiteGF(ind)%l.EQ.l) THEN
            EXIT
         ENDIF
         IF(ind.EQ.atoms%n_gf) CALL juDFT_error("Green's function element not found", hint="This is a bug in FLEUR, please report")
      ENDDO
   END SUBROUTINE

END MODULE m_hubbard1_setup