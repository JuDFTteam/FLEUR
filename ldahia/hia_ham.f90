MODULE m_hia_ham
   USE m_juDFT

   CONTAINS

   SUBROUTINE hia_setup(gOnsite,sym,atoms,sphhar, input,el,inDen,pot,mpi,results)

      !Sets up quantities needed for DFT+HIA (Has similar structure to u_setup.f90 but also sets up the stomic hamiltonians)

      USE m_umtx
      USE m_uj2f
      USE m_types
      USE m_vmmp
      USE m_constants
      USE m_fock_basis

      IMPLICIT NONE

      TYPE(t_greensf),INTENT(IN)      :: gOnsite
      TYPE(t_sym),INTENT(IN)          :: sym
      TYPE(t_results),INTENT(INOUT)   :: results
      TYPE(t_mpi),INTENT(IN)          :: mpi
      TYPE(t_input),INTENT(IN)        :: input
      TYPE(t_sphhar),INTENT(IN)       :: sphhar
      TYPE(t_atoms),INTENT(IN)        :: atoms
      TYPE(t_potden),INTENT(IN)       :: inDen
      TYPE(t_potden),INTENT(INOUT)    :: pot

      REAL,    INTENT(IN)           :: el(0:,:,:) !(0:atoms%lmaxd,ntype,jspd)

      INTEGER     i,i_hia,N,n_occ,i_gf,m,ind_h,N_basis,N_states
      INTEGER     itype,ispin,j,k,l,jspin,urec,i_u
      INTEGER     noded,nodeu,ios,lty(atoms%n_u)
      INTEGER     max_states, neig
      REAL        wronk, n_f, mu
      LOGICAL     n_exist
      CHARACTER*8 l_type*2,l_form*9

      REAL f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)

      REAL, ALLOCATABLE :: u(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: n_mmp(:,:,:,:)
      INTEGER, ALLOCATABLE :: basis(:)
      COMPLEX, ALLOCATABLE :: ev(:,:)
      COMPLEX, ALLOCATABLE :: eig(:)
      TYPE(t_mat) :: h_mat


      IF(ANY(gOnsite%gmmpMat(:,:,:,:,:,:).NE.CMPLX(0.0,0.0)).AND.atoms%n_hia.GT.0) THEN

         ! calculate slater integrals from u and j
         CALL uj2f(input%jspins,atoms%lda_hia(:),atoms%n_hia,f0,f2,f4,f6)

         ! set up e-e- interaction matrix
         ALLOCATE (u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                   -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_hia)))
         ALLOCATE (n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_hia),input%jspins))

         !If the slater integrals are calculated in cored they depend on spin
         !We average over them
         f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
         f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
         f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
         f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2

         CALL umtx(atoms%lda_hia(:),atoms%n_hia,f0(:,1),f2(:,1),f4(:,1),f6(:,1),u)
         !
         !Missing: matrix elements of the hamiltonian and soc 
         !

         !
         !Calculate the number of electrons in the correlated orbitals
         !
         CALL gOnsite%calc_mmpmat(atoms,sym,input%jspins,n_mmp)
         DO i_hia = 1, atoms%n_hia

            l = atoms%lda_hia(i_hia)%l
            n = atoms%lda_hia(i_hia)%atomType

            !Find the corresponding index of the onsite gf
            CALL gOnsite%index(l,n,i_gf)

            n_f = 0.0
            DO ispin = 1, input%jspins
               DO m = -l,l
                  n_f = n_f + n_mmp(m,m,i_gf,ispin)
               ENDDO
            ENDDO
            !
            !Set up and digaonalize the atomic hamiltonian for n-1, n and n+1 electrons in the correlated orbitals
            !n being the nearest integer value of n_f
            !
            n_occ = ANINT(n_f)

            DO N = n_occ-1, n_occ+1
               !TEMPORARY:
               mu = 0.0
               !
               !Find all 2*(2*l+1) bit numbers with N bits being 1
               !
               N_states = 2*(2*l+1) !number of one-particle states in the orbital
               CALL gen_fock_states(N,N_states,basis,N_basis)

               CALL h_mat%init(.true.,N_basis,N_basis)
               CALL hia_ham(l,N,u(:,:,:,:,i_hia),mu,basis(:),h_mat)
               !
               !Diagonalize the matrix
               !
               !At the moment we use the lapack routine (the goal should be to reuse the interface for diagonalization but we want to avoid creating the overlap matrix (unity in our case)) 
               !

               neig = 0
               eig = 0.0
               ev = 0.0
            ENDDO
            !
            !Calculate the interacting Green's function
            !






            !
            !Invert to obtain the self-energy
            !

            !
            !Add impurity to onsite Green's function
            !

            !
            !Calculate the density matrix and write it into potden%mmpmat
            !


         ENDDO
      ELSE 
         IF (mpi%irank.EQ.0) THEN
          WRITE (*,*) 'no onsite gf matrix found ... skipping LDA+Hubbard 1'
         ENDIF
      ENDIF

   END SUBROUTINE hia_setup


   SUBROUTINE add_self_energy(g0,e,g,mu_at,l)

      USE m_juDFT

      !This Subroutine calculates the (m,mp)-Matrix block for given site and spin of the interacting Green's function
      !This is done according to EQ(2) of the methods section in of Sci. Rep. 7, 2751 (2017)

      COMPLEX,         INTENT(IN)     :: g0(:,:) !non-interacting Green's function
      COMPLEX,         INTENT(IN)     :: e(:,:)  !self energy-correction from DFT+HIA
      COMPLEX,         INTENT(OUT)    :: g(:,:)  !interacting greens function
      REAL,            INTENT(IN)     :: mu_at   !this makes sure that the occupations do not deviate from g0 to g
      INTEGER,         INTENT(IN)     :: l      

      INTEGER matsize, i,j
      INTEGER              :: info
      INTEGER, ALLOCATABLE :: ipiv(:)
      COMPLEX, ALLOCATABLE :: work(:)
      COMPLEX, ALLOCATAbLE :: inv(:,:)

      matsize = 2*l+1

      ALLOCATE(work(matsize))
      ALLOCATE(inv(matsize,matsize))
      ALLOCATE(ipiv(matsize))

      !Array where to perform the operations
      inv(:,:) = g0(:,:)
      !
      !  invert the non-interacting greens function -Matrix
      !
      CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
      IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="add_self_energy")
      CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
      IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="add_self_energy")
      !
      !  add self energy correction
      !
      DO i = 1, matsize
         DO j = 1, matsize

            inv(i,j) = inv(i,j) + e(i,j)

            IF(i.EQ.j) inv(i,j) = inv(i,j) + mu_at

         ENDDO
      ENDDO
      !
      !  invert to obtain interacting greens function
      !
      CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
      IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="add_self_energy")
      CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
      IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="add_self_energy")

      g(:,:) = inv(:,:)

   END SUBROUTINE


   SUBROUTINE hia_ham(l,N,U,mu,basis,h_mat)

      USE m_types
      USE m_constants
      USE m_fock_basis

      !This is the naive approach to set up the Fock-Matrix for the atomic Hamiltonian

      INTEGER,                INTENT(IN)    :: l                                                     !l quantum number for the impurity
      INTEGER,                INTENT(IN)    :: N                                                     !number of electrons for which to set up the hamiltonian
      REAL,                   INTENT(IN)    :: mu                                                    !chemical potential
      REAL,                   INTENT(IN)    :: U(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&   
                                         -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)   !Full U-tensor for the site
      !REAL,                   INTENT(IN)    :: H(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)!Hamiltonian from DFT for the correlated orbitals (??)

      INTEGER,                INTENT(INOUT) :: basis(:)
      TYPE(t_mat),            INTENT(INOUT) :: h_mat

      INTEGER N_states, N_op, N_basis
      INTEGER dist
      INTEGER i,j,i_op
      INTEGER, ALLOCATABLE :: op(:,:)
      INTEGER :: m(4), s(4)


      !Idea: represent the fock states as binary numbers e.g 01001101001010 for a 6 electron state in the 4f orbitals
      !and use the corresponding integer to label them. This way every state is uniquely determined by the number 

      !Difficulty: keep track of what bit means what (done in m_fock_basis)
      !
      !At the moment we construct the full matrix and do not try to reduce space usage by using the sparseness
      !

      !
      !Add the chemical potential to the diagonal elements of the Fock-Matrix
      !
      DO i = 1, N_basis
         h_mat%data_r(i,i) = h_mat%data_r(i,i) + mu
      ENDDO
      !
      !Construct the matrix elements of the atomic hamiltonian
      !
      DO i = 1, N_basis
         DO j = i, N_basis
            !
            !Find the transitions contributing to the U-term (2 pairs of creation/annihilation and no spin change)
            !
            CALL find_transitions(basis(i),basis(j),N_states,2,.false.,N_op,op,dist)

            IF(N_op.NE.0) THEN
               !
               !Add the corresponding matrix elements of the U-tensor (be careful with the order)
               !
               DO i_op = 1, N_op
                  !
                  ! The operations are still encoded in bit indices (maybe change)
                  !
                  CALL bit_to_op(op(i_op,:),4,l,m(:),s(:))
                  !
                  ! The corresponding U-element is 1/2*U(m1, m2, m4, m3) 
                  !
                  h_mat%data_r(i,j) = h_mat%data_r(i,j) + 1/2. * U(m(1),m(2),m(4),m(3))
               ENDDO
            END IF
            !
            !Avoid second call to find_transitions with one pair and allowed spin change if the states are too different
            !
            IF(dist.LE.2) THEN
               CALL find_transitions(basis(i),basis(j),N_states,1,.true.,N_op,op,dist)

               IF(N_op.NE.0) THEN
                  !
                  !Add the corresponding matrix elements of the porjected hamitonian and other terms
                  !
                  DO i_op = 1, N_op
                     !
                     ! The operations are still encoded in bit indices (maybe change)
                     !
                     CALL bit_to_op(op(i_op,:),2,l,m(:),s(:))
                     !
                     ! 
                     !
                     !h_hia(i,j) = h_hia(i,j) + ???
                  ENDDO
               END IF
            END IF
         ENDDO
      ENDDO
   END SUBROUTINE hia_ham


   SUBROUTINE g_at_matrixelem(m,spin,l,nu,nuprime,b,bprime,N,Nprime)

      !This subroutine calculates the matrix elements needed for the interacting green's function 
      !of the atomic hamilatonian in the lehmann representation


   END SUBROUTINE

END MODULE m_hia_ham