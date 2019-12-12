MODULE m_hubbard1_setup

   USE m_juDFT

   IMPLICIT NONE

   LOGICAL, PARAMETER :: l_setupdebug = .FALSE.  !Enable/Disable Debug outputs like dependency of occupation on chemical potential shift 
   CHARACTER(len=30), PARAMETER :: main_folder = "Hubbard1"



   CONTAINS

   SUBROUTINE hubbard1_setup(atoms,input,sym,mpi,noco,pot,gdft,hub1,results,den)

      USE m_types
      USE m_constants
      USE m_uj2f
      USE m_mudc
      USE m_denmat_dist
      USE m_gfcalc
      USE m_hubbard1_io
      USE m_add_selfen
#ifdef CPP_EDSOLVER
      USE EDsolver, only: EDsolver_from_cfg
#endif
#ifdef CPP_MPI
      INCLUDE "mpif.h"
#endif


      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_mpi),      INTENT(IN)     :: mpi
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_potden),   INTENT(IN)     :: pot
      TYPE(t_greensf),  INTENT(IN)     :: gdft !green's function calculated from the Kohn-Sham system
      TYPE(t_hub1ham),  INTENT(INOUT)  :: hub1
      TYPE(t_results),  INTENT(INOUT)  :: results
      TYPE(t_potden),   INTENT(INOUT)  :: den

#ifdef CPP_MPI
      EXTERNAL MPI_BCAST
#endif

      !-- Local Scalars
      INTEGER :: i_hia,nType,l,n_occ,ispin,m,iz,k,j,i_exc,i,jspin,ipm
      INTEGER :: io_error,ierr
      INTEGER :: indStart,indEnd
      REAL    :: mu_dc,exc
      LOGICAL :: l_selfenexist,l_exist,l_linkedsolver,l_ccfexist,l_bathexist,occ_err

      CHARACTER(len=300) :: cwd,path,folder,xPath
      CHARACTER(len=8)   :: l_type*2,l_form*9

      !-- Local Types
      TYPE(t_greensf)    :: gu

      !-- Local Arrays
      REAL    :: f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL    :: f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      REAL    :: n_l(atoms%n_hia,input%jspins)
      COMPLEX :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,3)
      COMPLEX :: e(gdft%nz)
      COMPLEX, ALLOCATABLE :: selfen(:,:,:,:,:)

      !To avoid confusing structure with pre-processor-switches later on
#ifdef CPP_EDSOLVER
      l_linkedsolver = .TRUE.
#else
      l_linkedsolver = .FALSE.
#endif

      e = 0.0

      !Positions of the DFT+HIA elements in all DFT+U related arrays
      indStart = atoms%n_u+1
      indEnd   = atoms%n_u+atoms%n_hia


      IF(mpi%irank.EQ.0) THEN
         !The onsite green's function was calculated but the solver
         !was not yet run
         !--> write out the configuration for the hubbard 1 solver
         IF(.NOT.ANY(gdft%gmmpMat(:,:,:,:,:,:).NE.0.0)) CALL juDFT_error("Hubbard-1 has no DFT greensf available",calledby="hubbard1_setup")
         CALL gu%init(input,lmaxU_const,atoms,noco,nz_in=gdft%nz, e_in=gdft%e,de_in=gdft%de,matsub_in=gdft%nmatsub)
         ALLOCATE(selfen(2*(2*lmaxU_const+1),2*(2*lmaxU_const+1),gdft%nz,2,atoms%n_hia))
         selfen = 0.0


         !Get the working directory
         CALL get_environment_variable('PWD',cwd)
         path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(main_folder))
         CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)))
         !Remove everything from the last iteration (Good Idea??)
         CALL SYSTEM('rm -rf ' // TRIM(ADJUSTL(path)) // "/*")

         ! calculate slater integrals from u and j
         CALL uj2f(input%jspins,atoms%lda_u(indStart:indEnd),atoms%n_hia,f0,f2,f4,f6)

         DO ispin = 1, 1 ! input%jspins
            f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
            f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
            f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
            f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
         END DO


         DO i_hia = 1, atoms%n_hia
            nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
            l = atoms%lda_u(atoms%n_u+i_hia)%l

            CALL hubbard1_path(atoms,i_hia,folder)
            WRITE(xPath,*) TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(xPath)))

            !For the first iteration we can fix the occupation and magnetic moments in the inp.xml file
            IF(hub1%iter.EQ.1.AND.ALL(den%mmpMat(:,:,indStart:indEnd,:).EQ.0.0)) THEN
               n_l(i_hia,:) = hub1%init_occ(i_hia)/input%jspins
               DO i_exc = 1, hub1%n_exc_given(i_hia)
                  hub1%mag_mom(i_hia,i_exc) = hub1%init_mom(i_hia,i_exc)
               ENDDO
            ELSE
               !calculate the occupation of the correlated shell
               CALL occmtx(gdft,l,nType,atoms,sym,input,mmpMat(:,:,i_hia,:),occ_err)
               n_l(i_hia,:) = 0.0
               DO ispin = 1, input%jspins
                  DO m = -l, l
                     n_l(i_hia,ispin) = n_l(i_hia,ispin) + REAL(mmpMat(m,m,i_hia,ispin))
                  ENDDO
               ENDDO
            ENDIF


            !Initial Information (We are already on irank 0)
            WRITE(6,*)
            WRITE(6,9010) nType
            WRITE(6,"(A)") "Everything related to the solver (e.g. mu_dc) is given in eV"
            WRITE(6,*)
            WRITE(6,"(A)") "Occupation from DFT-Green's function:"
            WRITE(6,9020) 'spin-up','spin-dn'
            WRITE(6,9030) n_l(i_hia,:)

            !--------------------------------------------------------------------------
            ! Calculate the chemical potential for the solver
            ! This is equal to the double-counting correction used in DFT+U
            !--------------------------------------------------------------------------
            ! V_FLL = U (n - 1/2) - J (n - 1) / 2
            ! V_AMF = U n/2 + 2l/[2(2l+1)] (U-J) n
            !--------------------------------------------------------------------------
            CALL mudc(atoms%lda_u(atoms%n_u+i_hia),n_l(i_hia,:),input%jspins,mu_dc)

            !Check wether the hubbard 1 solver was run:(old version)
            INQUIRE(file=TRIM(ADJUSTL(xPath)) // "se.atom",exist=l_selfenexist)

            IF(mpi%irank.EQ.0.AND..NOT.l_selfenexist) THEN
               !Nearest Integer occupation
               n_occ = ANINT(SUM(n_l(i_hia,:)))
               IF(l_linkedsolver) THEN
                  !-------------------------------------------------------
                  ! Check for additional input files
                  !-------------------------------------------------------
                  !Is a crystal field matrix present in the work directory
                  INQUIRE(file=TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(cfg_file_ccf)),exist=l_ccfexist)
                  IF(l_ccfexist) CALL read_ccfmat(TRIM(ADJUSTL(cwd)),hub1%ccfmat(i_hia,-l:l,-l:l),l)
                  !Is a bath parameter file present
                  INQUIRE(file=TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(cfg_file_bath)),exist=l_bathexist)
                  !Copy the bath file to the Hubbard 1 solver if its present
                  IF(l_bathexist) CALL SYSTEM('cp ' // TRIM(ADJUSTL(cfg_file_bath)) // ' ' // TRIM(ADJUSTL(xPath)))
                  !-------------------------------------------------------
                  ! Write the main config files
                  !-------------------------------------------------------
                  CALL hubbard1_input(xPath,i_hia,l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),hub1,mu_dc,n_occ,&
                                     l_bathexist,hub1%iter==1.AND.ALL(pot%mmpmat(:,:,indStart:indEnd,:).EQ.0.0),.true.)
               ELSE
                  !If no Solver is linked we assume that the old solver is used and we write out some additional files
                  !Crystal field matrix (old version)
                  IF(hub1%ccf(i_hia).NE.0.0) THEN
                     CALL write_ccfmat(xPath,hub1%ccfmat(i_hia,-l:l,-l:l),l)
                  ENDIF
                  !Energy contour (old version)
                  IF(gdft%mode.NE.3) THEN
                     OPEN(unit=1337,file=TRIM(ADJUSTL(xPath)) // "contour.dat",status="replace",action="write")
                     DO iz = 1, gdft%nz
                        WRITE(1337,"(2f14.8)") REAL(gdft%e(iz))*hartree_to_ev_const, AIMAG(gdft%e(iz))*hartree_to_ev_const
                     ENDDO
                     CLOSE(unit= 1337)
                  ENDIF
                  !-------------------------------------------------------
                  ! Write the main config files (old version)
                  !-------------------------------------------------------
                  CALL hubbard1_input(xPath,i_hia,l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),hub1,mu_dc,n_occ,.false.,.false.,.false.)
               ENDIF
            ENDIF
         ENDDO

         !------------------------------------------------------------
         ! This loop runs the solver if it is available
         ! If not the program terminates here
         !------------------------------------------------------------
         DO i_hia = 1, atoms%n_hia
            nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
            l = atoms%lda_u(atoms%n_u+i_hia)%l

            CALL hubbard1_path(atoms,i_hia,folder)
            WRITE(xPath,*) TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
            INQUIRE(file=TRIM(ADJUSTL(xPath)) // "se.atom",exist=l_selfenexist)
            IF(mpi%irank.EQ.0) THEN
               IF(.NOT.l_selfenexist) THEN
                  IF(.NOT.l_linkedsolver) CALL juDFT_END("Hubbard1 input has been written into Hubbard1/ (No Solver linked)",mpi%irank)
                  CALL timestart("Hubbard 1: EDsolver")
                  !We have to change into the Hubbard1 directory so that the solver routines can read the config
                  CALL CHDIR(TRIM(ADJUSTL(xPath)))
#ifdef CPP_EDSOLVER
                  e = gdft%e*hartree_to_ev_const
                  CALL EDsolver_from_cfg(2*(2*l+1),gdft%nz,e,selfen(:,:,:,1,i_hia),1)
                  !---------------------------------------------------
                  ! Calculate selfenergy on lower contour explicitly
                  ! Mainly out of paranoia :D
                  ! No rediagonalization (last argument switches this)
                  !---------------------------------------------------
                  e = conjg(gdft%e)*hartree_to_ev_const
                  CALL EDsolver_from_cfg(2*(2*l+1),gdft%nz,e,selfen(:,:,:,2,i_hia),0)
#endif
                  CALL CHDIR(TRIM(ADJUSTL(cwd)))
                  CALL timestop("Hubbard 1: EDsolver")
               ELSE
                  !If there is no linked solver library we read in the selfenergy here
                  CALL read_selfen(xPath,selfen(:,:,:,1,i_hia),gdft%nz,2*(2*l+1),.false.)
               ENDIF
            ENDIF
         ENDDO

         IF(l_selfenexist.OR.l_linkedsolver) THEN
            !-------------------------------------------
            ! Postprocess selfenergy
            !-------------------------------------------
            DO i_hia = 1, atoms%n_hia
               nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
               l = atoms%lda_u(atoms%n_u+i_hia)%l
               DO ipm = 1, 2
                  DO iz = 1, gdft%nz
                     !---------------------------------------------
                     ! Convert the selfenergy to hartree
                     !---------------------------------------------
                     selfen(:,:,iz,ipm,i_hia) = selfen(:,:,iz,ipm,i_hia)/hartree_to_ev_const
                     !---------------------------------------------
                     ! The order of spins is reversed in the Solver
                     !---------------------------------------------
                     CALL swapSpin(selfen(:,:,iz,ipm,i_hia),2*l+1)
                     !---------------------------------------------------------------------
                     ! The DFT green's function also includes the previous DFT+U correction
                     ! This is removed by substracting it from the selfenergy
                     !---------------------------------------------------------------------
                     CALL removeU(selfen(:,:,iz,ipm,i_hia),l,input%jspins,noco%l_mtNocoPot,pot%mmpMat(:,:,atoms%n_u+i_hia,:))
                  ENDDO
               ENDDO
            ENDDO
            !----------------------------------------------------------------------
            ! Solution of the Dyson Equation
            !----------------------------------------------------------------------
            ! G(z)^(-1) = G_0(z)^(-1) - mu - Sigma(z)
            !----------------------------------------------------------------------
            ! Sigma(z) is the self-energy from the impurity solver
            ! We introduce an additional chemical potential mu, which is determined
            ! so that the occupation of the correlated orbital does not change
            !----------------------------------------------------------------------
            CALL timestart("Hubbard 1: Add Selfenenergy")
            CALL add_selfen(gdft,selfen,atoms,input,sym,noco,hub1,results%ef,n_l,gu,mmpMat)
            CALL timestop("Hubbard 1: Add Selfenenergy")
            IF(l_setupdebug) THEN
               DO i_hia = 1, atoms%n_hia
                  nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
                  l = atoms%lda_u(atoms%n_u+i_hia)%l
                  CALL gfDOS(gdft,l,nType,800+i_hia+hub1%iter,atoms,input,results%ef)
                  CALL gfDOS(gu,l,nType,900+i_hia+hub1%iter,atoms,input,results%ef)
                  CALL writeSelfenElement(selfen(:,:,:,1,i_hia),gdft%e,results%ef,gdft%nz,2*l+1)
               ENDDO
            ENDIF
            !----------------------------------------------------------------------
            ! Calculate the distance and update the density matrix
            !----------------------------------------------------------------------
            DO i_hia = 1, atoms%n_hia
               CALL n_mmp_dist(den%mmpMat(:,:,atoms%n_u+i_hia,:),mmpMat(:,:,i_hia,:),input,results)
               DO ispin = 1, MERGE(3,input%jspins,noco%l_mperp)
                  den%mmpMat(:,:,atoms%n_u+i_hia,ispin) = mmpMat(:,:,i_hia,ispin)
               ENDDO
            ENDDO
         ENDIF
         DEALLOCATE(selfen)
      ELSE
         !If we are on a different mpi%irank and no solver is linked we need to call juDFT_end here if the solver was not run
         !kind of a weird workaround (replace with something better)
         IF(.NOT.l_linkedsolver) THEN
            CALL get_environment_variable('PWD',cwd)
            DO i_hia = atoms%n_u+1, atoms%n_u+atoms%n_hia
               CALL hubbard1_path(atoms,i_hia,folder)
               WRITE(xPath,*) TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
               INQUIRE(file=TRIM(ADJUSTL(xPath)) // "se.atom",exist=l_selfenexist)
               IF(.NOT.l_selfenexist) CALL juDFT_END("Hubbard1 input has been written into Hubbard1/ (No Solver linked)",mpi%irank)
            ENDDO
            !If we are here the solver was run and we go to MPI_BCAST
         ENDIF
      ENDIF
#ifdef CPP_MPI
      !Broadcast the density matrix here
      CALL MPI_BCAST(den%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,indStart:indEnd,:),&
                     49*atoms%n_hia*MERGE(3,input%jspins,noco%l_mperp),MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
#endif

      !FORMAT Statements:
9010  FORMAT("Setup for Hubbard 1 solver for atom ", I3, ": ")
9020  FORMAT(TR8,A7,TR3,A7)
9030  FORMAT(TR7,f8.4,TR2,f8.4)

   END SUBROUTINE hubbard1_setup

   SUBROUTINE swapSpin(mat,ns)

      IMPLICIT NONE

      COMPLEX,       INTENT(INOUT) :: mat(:,:)
      INTEGER,       INTENT(IN)    :: ns

      COMPLEX tmp(2*ns,2*ns),tmp_off(ns,ns)
      INTEGER i,j

      !Transformation matrix
      tmp = 0.0
      tmp_off = 0.0
      DO i = 1, ns
         tmp(i,ns+i) = 1.0
         tmp(ns+i,i) = 1.0
         tmp_off(i,ns-i+1) = 1.0
      ENDDO
      !WRITE(*,*) "BEFORE"
      !WRITE(*,"(14f8.5)") REAL(mat)
      mat = matmul(mat,tmp)
      mat = matmul(tmp,mat)

      !mat(1:ns,ns+1:2*ns) = matmul(mat(1:ns,ns+1:2*ns),tmp_off)
      !mat(1:ns,ns+1:2*ns) = transpose(matmul(tmp_off,mat(1:ns,ns+1:2*ns)))

      !mat(ns+1:2*ns,1:ns) = matmul(mat(ns+1:2*ns,1:ns),tmp_off)
      !mat(ns+1:2*ns,1:ns) = transpose(matmul(tmp_off,mat(ns+1:2*ns,1:ns)))

      !WRITE(*,*) "AFTER"
      !WRITE(*,"(14f8.5)") REAL(mat)

   END SUBROUTINE swapSpin


   SUBROUTINE removeU(mat,l,jspins,l_vmperp,vmmp)

      USE m_constants

      IMPLICIT NONE
      COMPLEX,       INTENT(INOUT) :: mat(:,:)
      INTEGER,       INTENT(IN)    :: l
      INTEGER,       INTENT(IN)    :: jspins
      LOGICAL,       INTENT(IN)    :: l_vmperp
      COMPLEX,       INTENT(IN)    :: vmmp(-lmaxU_const:,-lmaxU_const:,:)

      INTEGER ns,i,j,m,mp,ispin

      ns = 2*l+1

      DO i = 1, ns
         DO j = 1, ns
            m  = i-1-l
            mp = j-1-l
            DO ispin = 1, MERGE(3,jspins,l_vmperp)
               IF(ispin < 3) THEN
                  mat(i+(ispin-1)*ns,j+(ispin-1)*ns) = mat(i+(ispin-1)*ns,j+(ispin-1)*ns) - REAL(vmmp(m,mp,ispin))/(3.0-jspins)
                  IF(jspins.EQ.1) mat(i+ns,j+ns) = mat(i+ns,j+ns) - REAL(vmmp(-m,-mp,ispin))/(3.0-jspins)
               ELSE
                  !----------------------------------------------------------------------------
                  ! The offdiagonal elements only have to be removed if they are actually added
                  ! to the hamiltonian (so noco%l_mperp and noco%l_mtNocoPot)
                  !----------------------------------------------------------------------------
                  mat(i+ns,j) = mat(i+ns,j) - vmmp(m,mp,ispin)
                  mat(i,j+ns) = mat(i,j+ns) - conjg(vmmp(mp,m,ispin))
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE removeU


   SUBROUTINE hubbard1_path(atoms,i_hia,xPath)

      !Defines the folder structure
      ! The Solver is run in the subdirectories
      ! Hubbard1/ if only one Hubbard1 prodcedure is run
      ! Hubbard1/atom_label_l if there are more

      USE m_types

      IMPLICIT NONE

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      INTEGER,             INTENT(IN)  :: i_hia
      CHARACTER(len=300),  INTENT(OUT) :: xPath

      CHARACTER(len=300) :: folder,fmt
      CHARACTER(len=1)   :: l_name(0:3)
      INTEGER nType,l

      l_name(0:3) = (/"s","p","d","f"/)

      nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
      l = atoms%lda_u(atoms%n_u+i_hia)%l
      xPath = TRIM(ADJUSTL(main_folder))
      IF(atoms%n_hia>1) THEN
         WRITE(fmt,'("(A",I2.2,",A1,A1,A1)")') LEN(TRIM(ADJUSTL(atoms%label(nType))))
         WRITE(folder,fmt) TRIM(ADJUSTL(atoms%label(nType))),"_",l_name(l),"/"
      ELSE
         folder=""
      ENDIF
      xPath = TRIM(ADJUSTL(xPath)) // "/" // folder

   END SUBROUTINE hubbard1_path

   SUBROUTINE writeSelfenElement(selfen,e,ef,nz,ns)

      USE m_constants

      IMPLICIT NONE

      INTEGER,       INTENT(IN)  :: nz,ns
      REAL,          INTENT(IN)  :: ef
      COMPLEX,       INTENT(IN)  :: selfen(:,:,:)
      COMPLEX,       INTENT(IN)  :: e(:)

      INTEGER iz,i,io_error
      CHARACTER(len=300) file
      COMPLEX up, down

      OPEN(unit=3456,file="selfen",status="replace",iostat = io_error)

      DO iz = 1, nz
         up = 0.0
         down = 0.0
         DO i = 1, ns
            up = up + selfen(i,i,iz)
            down = down + selfen(i+ns,i+ns,iz)
         ENDDO
         WRITE(3456,"(5f14.8)") REAL(e(iz))*hartree_to_ev_const, -AIMAG(up), -AIMAG(down), REAL(up), REAL(down)
      ENDDO

      CLOSE(unit=3456)


   END SUBROUTINE writeSelfenElement

END MODULE m_hubbard1_setup
