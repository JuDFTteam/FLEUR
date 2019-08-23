MODULE m_hubbard1_setup

   USE m_juDFT

   IMPLICIT NONE

   LOGICAL, PARAMETER :: l_debug = .TRUE.  !Enable/Disable Debug outputs like dependency of occupation on chemical potential shift 
   CHARACTER(len=30), PARAMETER :: main_folder = "Hubbard1"



   CONTAINS

   SUBROUTINE hubbard1_setup(atoms,hub1,sym,mpi,noco,input,usdus,den,pot,gdft,results)

      USE m_types
      USE m_hubbard1_io
      USE m_uj2f
      USE m_constants
      USE m_gfcalc
      USE m_umtx
      USE m_vmmp
      USE m_add_selfen
      USE m_mudc
      USE m_occmtx
      USE m_denmat_dist
#ifdef CPP_EDSOLVER
      USE EDsolver, only: EDsolver_from_cfg
#endif
#ifdef CPP_MPI
      INCLUDE "mpif.h"
#endif


      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_hub1ham),  INTENT(INOUT)  :: hub1
      TYPE(t_sym),      INTENT(IN)     :: sym
      TYPE(t_mpi),      INTENT(IN)     :: mpi
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_usdus),    INTENT(IN)     :: usdus
      TYPE(t_potden),   INTENT(INOUT)  :: den
      TYPE(t_potden),   INTENT(INOUT)  :: pot
      TYPE(t_results),  INTENT(INOUT)  :: results
      TYPE(t_greensf),  INTENT(IN)     :: gdft !green's function calculated from the Kohn-Sham system
      
#ifdef CPP_MPI
      EXTERNAL MPI_BCAST
#endif

      INTEGER i_hia,nType,l,n_occ,ispin,m,iz,k,j,i_exc,i
      INTEGER io_error,ierr
      INTEGER indStart,indEnd
      REAL    mu_dc,e_lda_hia,exc

      CHARACTER(len=300) :: cwd,path,folder,xPath
      CHARACTER(len=8)   :: l_type*2,l_form*9

      TYPE(t_greensf)    :: gu 

      REAL     f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL     f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      REAL     u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
               -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia)

      COMPLEX  mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,3)
      COMPLEX  selfen(atoms%n_hia,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1),gdft%nz)
      REAL     n_l(atoms%n_hia,input%jspins)
      LOGICAL  l_selfenexist,l_exist,l_linkedsolver,l_ccfexist,l_bathexist

      !Positions of the DFT+HIA elements in all DFT+U related arrays
      indStart = atoms%n_u+1
      indEnd   = atoms%n_u+atoms%n_hia 


      !TODO: We don't need to calculate the green's function in every iteration 
      !(what is an appropriate cutoff for the distance under which we calculate the greens function)

      IF(ANY(den%mmpMat(:,:,indStart:indEnd,:).NE.0.0).OR.hub1%l_runthisiter) THEN
         !Get the slater integrals from the U and J parameters
         CALL uj2f(input%jspins,atoms%lda_u(indStart:indEnd),atoms%n_hia,f0,f2,f4,f6)
         f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
         f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
         f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
         f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
         CALL umtx(atoms%lda_u(indStart:indEnd),atoms%n_hia,f0(:,1),f2(:,1),f4(:,1),f6(:,1),u)

         CALL v_mmp(sym,atoms,atoms%lda_u(indStart:indEnd),atoms%n_hia,input%jspins,.FALSE.,den%mmpMat(:,:,indStart:indEnd,:),&
         u,f0,f2,pot%mmpMat(:,:,indStart:indEnd,:),e_lda_hia)

         IF(hub1%l_runthisiter.AND.(ANY(gdft%gmmpMat(:,:,:,:,:,:).NE.0.0)).AND.mpi%irank.EQ.0) THEN 
            !The onsite green's function was calculated but the solver 
            !was not yet run
            !--> write out the configuration for the hubbard 1 solver 
            CALL gu%init(input,lmaxU_const,atoms,noco,nz_in=gdft%nz, e_in=gdft%e,de_in=gdft%de,matsub_in=gdft%nmatsub)

            !Get the working directory
            CALL get_environment_variable('PWD',cwd)
            path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(main_folder))
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)))
            !Remove everything from the last iteration (Good Idea??)
            CALL SYSTEM('rm -rf ' // TRIM(ADJUSTL(path)) // "/*")

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
                  CALL occmtx(gdft,l,nType,atoms,sym,input,mmpMat(:,:,i_hia,:))
                  n_l(i_hia,:) = 0.0
                  DO ispin = 1, input%jspins
                     DO m = -l, l
                        n_l(i_hia,ispin) = n_l(i_hia,ispin) + mmpMat(m,m,i_hia,ispin)
                     ENDDO
                  ENDDO
               ENDIF

               !Initial Information
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
               CALL mudc(atoms%lda_u(atoms%n_u+i_hia),n_l(i_hia,:),mu_dc,input%jspins)

               !Check wether the hubbard 1 solver was run:(old version)
               INQUIRE(file=TRIM(ADJUSTL(xPath)) // "se.atom",exist=l_selfenexist)

               IF(mpi%irank.EQ.0.AND..NOT.l_selfenexist) THEN
                  !Nearest Integer occupation
                  n_occ = ANINT(SUM(n_l(i_hia,:)))
#ifdef CPP_EDSOLVER
                  !-------------------------------------------------------
                  ! Check for additional input files
                  !-------------------------------------------------------
                  !Is a crystal field matrix present in the work directory
                  WRITE(*,*) TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(cfg_file_ccf))
                  INQUIRE(file=TRIM(ADJUSTL(cwd)) // TRIM(ADJUSTL(cfg_file_ccf)),exist=l_ccfexist)
                  IF(l_ccfexist) CALL read_ccfmat(TRIM(ADJUSTL(cwd)),hub1%ccfmat(i_hia,-l:l,-l:l),l)
                  !Is a bath parameter file present 
                  INQUIRE(file=TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(cfg_file_bath)),exist=l_bathexist)
                  !Copy the bath file to the Hubbard 1 solver if its present
                  IF(l_bathexist) CALL SYSTEM('cp ' // TRIM(ADJUSTL(cfg_file_bath)) // ' ' // TRIM(ADJUSTL(xPath)))
                  !-------------------------------------------------------
                  ! Write the main config files
                  !-------------------------------------------------------
                  CALL hubbard1_input(xPath,i_hia,l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),hub1,mu_dc,n_occ,l_bathexist,.true.)
#else
                  CALL hubbard1_input(xPath,i_hia,l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),hub1,mu_dc,n_occ,.false.,.false.)

                  !If no Solver is linked we assume that the old solver is used and we write out some additional files
                  !Crystal field matrix (old version)
                  IF(hub1%ccf(i_hia).NE.0.0) THEN
                     CALL write_ccfmat(xPath,hub1%ccfmat(i_hia,-l:l,-l:l),l)
                  ENDIF
                  !Energy contour (old version)
                  IF(gdft%mode.NE.1) THEN
                     OPEN(unit=1337,file=TRIM(ADJUSTL(xPath)) // "contour.dat",status="replace",action="write")
                     DO iz = 1, gdft%nz
                        WRITE(1337,"(2f14.8)") REAL(gdft%e(iz))*hartree_to_ev_const, AIMAG(gdft%e(iz))*hartree_to_ev_const
                     ENDDO
                     CLOSE(unit= 1337)
                  ENDIF   
#endif
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
#ifdef CPP_EDSOLVER
                     l_linkedsolver=.TRUE.
                     CALL timestart("Hubbard 1: EDsolver")
                     !We have to change into the Hubbard1 directory so that the solver routines can read the config
                     CALL CHDIR(TRIM(ADJUSTL(xPath)))                  
                     CALL EDsolver_from_cfg(2*(2*l+1),gdft%nz,gdft%e(1:gdft%nz)*hartree_to_ev_const,selfen(i_hia,:,:,1:gdft%nz),1)
                     !The solver is given everything in eV by default, so we need to convert back to htr
                     selfen(i_hia,:,:,:) = selfen(i_hia,:,:,:)/hartree_to_ev_const
                     CALL CHDIR(TRIM(ADJUSTL(cwd)))
                     CALL timestop("Hubbard 1: EDsolver")
#else 
                     CALL juDFT_END("Hubbard1 input has been written into Hubbard1/ (No Solver linked)",mpi%irank)
#endif
                  ELSE               
                     !If there is no linked solver library we read in the selfenergy here
                     CALL read_selfen(xPath,selfen(i_hia,1:2*(2*l+1),1:2*(2*l+1),1:gdft%nz),gdft%nz,2*(2*l+1),.false.)
                  ENDIF
               ENDIF
            ENDDO

            IF(l_selfenexist.OR.l_linkedsolver) THEN
               DO i_hia = 1, atoms%n_hia
                  nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
                  l = atoms%lda_u(atoms%n_u+i_hia)%l
                  DO iz = 1, gdft%nz
                     CALL swapSpin(selfen(i_hia,:,:,iz),2*l+1)
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
               CALL add_selfen(gdft,gu,selfen,atoms,noco,hub1,sym,input,results%ef,n_l,mu_dc/hartree_to_ev_const,&
                              pot%mmpMat(:,:,indStart:indEnd,:),mmpMat)
               CALL timestop("Hubbard 1: Add Selfenenergy")
               !----------------------------------------------------------------------
               ! Calculate DFT+U potential correction
               !----------------------------------------------------------------------
               CALL v_mmp(sym,atoms,atoms%lda_u(indStart:indEnd),atoms%n_hia,input%jspins,.FALSE.,mmpMat,&
                     u,f0,f2,pot%mmpMat(:,:,indStart:indEnd,:),results%e_ldau)
               !----------------------------------------------------------------------
               ! Calculate the distance and update the density matrix 
               !----------------------------------------------------------------------
               CALL n_mmp_dist(den%mmpMat(:,:,indStart:indEnd,:),mmpMat,atoms%n_hia,results,input%jspins)
               den%mmpMat(:,:,indStart:indEnd,:) = mmpMat(:,:,:,1:input%jspins) !For now LDA+U in FLEUR ignores spin offdiagonal elements
            ENDIF
         ELSE 
            !The solver does not need to be run so we just add the current energy correction from LDA+HIA 
            results%e_ldau = results%e_ldau + e_lda_hia 
         ENDIF
         !Write out the density matrix and potential matrix (compare u_setup.f90)
         IF (mpi%irank.EQ.0) THEN
            DO ispin = 1,input%jspins
               WRITE (6,'(a7,i3)') 'spin #',ispin
               DO i_hia = atoms%n_u+1, atoms%n_u+atoms%n_hia
                  nType = atoms%lda_u(i_hia)%atomType
                  l = atoms%lda_u(i_hia)%l
                  WRITE (l_type,'(i2)') 2*(2*l+1)
                  l_form = '('//l_type//'f12.7)'
                  WRITE (6,'(a20,i3)') 'n-matrix for atom # ',nType
                  WRITE (6,l_form) ((den%mmpMat(k,j,i_hia,ispin),k=-l,l),j=-l,l)
                  WRITE (6,'(a20,i3)') 'V-matrix for atom # ',nType
                  IF (atoms%lda_u(i_hia)%l_amf) THEN
                     WRITE (6,*) 'using the around-mean-field limit '
                  ELSE
                     WRITE (6,*) 'using the atomic limit of LDA+U '
                  ENDIF
                  WRITE (6,l_form) ((pot%mmpMat(k,j,i_hia,ispin),k=-l,l),j=-l,l)
               END DO
            END DO
            WRITE (6,*) results%e_ldau
         ENDIF
      ELSE IF(mpi%irank.NE.0) THEN
         pot%mmpMat(:,:,atoms%n_u+1:atoms%n_hia+atoms%n_u,:) = CMPLX(0.0,0.0)
         !If we are on a different mpi%irank and no solver is linked we need to call juDFT_end here if the solver was not run
         !kind of a weird workaround (replace with something better)
#ifdef CPP_EDSOLVER 
         !Do nothing and go to the MPI_BCAST
#else 
         CALL get_environment_variable('PWD',cwd)
         DO i_hia = atoms%n_u+1, atoms%n_u+atoms%n_hia
            CALL hubbard1_path(atoms,i_hia,folder)
            WRITE(xPath,*) TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder)) 
            INQUIRE(file=TRIM(ADJUSTL(xPath)) // "se.atom",exist=l_selfenexist)
            IF(.NOT.l_selfenexist) CALL juDFT_END("Hubbard1 input has been written into Hubbard1/ (No Solver linked)",mpi%irank)
         ENDDO
         !If we are here the solver was run and we go to MPI_BCAST
#endif
      ELSE
         !occupation matrix is zero and LDA+Hubbard 1 shouldn't be run yet
         !There is nothing to be done yet just set the potential correction to 0
         pot%mmpMat(:,:,atoms%n_u+1:atoms%n_hia+atoms%n_u,:) = CMPLX(0.0,0.0)
         results%e_ldau = 0.0
      ENDIF
#ifdef CPP_MPI
      !Broadcast both the potential and the density matrix here
      CALL MPI_BCAST(pot%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,indStart:indEnd,:),&
                     49*atoms%n_hia*input%jspins,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
      CALL MPI_BCAST(den%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,indStart:indEnd,:),&
                     49*atoms%n_hia*input%jspins,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
#endif 
      
      !FORMAT Statements:
9010  FORMAT("Setup for Hubbard 1 solver for atom ", I3, ": ")
9020  FORMAT(TR8,A7,TR3,A7)
9030  FORMAT(TR7,f8.4,TR2,f8.4)

   END SUBROUTINE hubbard1_setup

   SUBROUTINE swapSpin(mat,ns)

      IMPLICIT NONE

      COMPLEX,       INTENT(INOUT) :: mat(2*ns,2*ns)
      INTEGER,       INTENT(IN)    :: ns

      COMPLEX tmp(ns,ns)

      !Spin-diagonal
      tmp = mat(1:ns,1:ns)
      mat(1:ns,1:ns) = mat(ns+1:2*ns,ns+1:2*ns)
      mat(ns+1:2*ns,ns+1:2*ns) = tmp
      !Spin-offdiagonal
      tmp = mat(ns+1:2*ns,1:ns)
      mat(ns+1:2*ns,1:ns) = transpose(conjg(mat(1:ns,ns+1:2*ns)))
      mat(1:ns,ns+1:2*ns) = transpose(conjg(tmp))

   END SUBROUTINE


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
         WRITE(fmt,'("(A",I2.2,",A1,A1,A1")') LEN(TRIM(ADJUSTL(atoms%label(nType))))
         WRITE(folder,fmt) atoms%label(nType),"_",l_name(l),"/"
      ELSE
         folder=""
      ENDIF
      xPath = TRIM(ADJUSTL(xPath)) // "/" // folder

   END SUBROUTINE hubbard1_path

   SUBROUTINE writeSelfenElement(selfen,e,nz,matsize,i)

      IMPLICIT NONE

      INTEGER,       INTENT(IN)  :: nz,matsize,i
      COMPLEX,       INTENT(IN)  :: selfen(matsize,matsize,nz)
      COMPLEX,       INTENT(IN)  :: e(nz)

      INTEGER iz
      CHARACTER(len=300) file

3456  FORMAT("selfen.",I3)
      WRITE(file,3456) i

      OPEN(unit=3456,file=file,status="replace")

      DO iz = 1, nz
         WRITE(3456,"(3f14.8)") REAL(e(iz)), REAL(selfen(i,i,iz)), AIMAG(selfen(i,i,iz))
      ENDDO

      CLOSE(unit=3456)


   END SUBROUTINE writeSelfenElement

END MODULE m_hubbard1_setup