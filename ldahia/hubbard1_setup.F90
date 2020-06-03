MODULE m_hubbard1_setup

   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_uj2f
   USE m_mudc
   USE m_hubbard1Distance
   USE m_occmtx
   USE m_hubbard1_io
   USE m_types_selfen
   USE m_add_selfen
   USE m_mpi_bc_tool
   USE m_greensf_io
#ifdef CPP_EDSOLVER
   USE EDsolver, only: EDsolver_from_cfg
#endif

   IMPLICIT NONE

#ifdef CPP_MPI
   INCLUDE 'mpif.h'
#endif
#include"cpp_double.h"

   CHARACTER(len=30), PARAMETER :: hubbard1CalcFolder = "Hubbard1"
   CHARACTER(len=30), PARAMETER :: hubbard1Outfile    = "out"

   CONTAINS

   SUBROUTINE hubbard1_setup(atoms,gfinp,hub1inp,input,mpi,noco,pot,gdft,hub1data,results,den)

      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_gfinp),    INTENT(IN)     :: gfinp
      TYPE(t_hub1inp),  INTENT(IN)     :: hub1inp
      TYPE(t_input),    INTENT(IN)     :: input
      TYPE(t_mpi),      INTENT(IN)     :: mpi
      TYPE(t_noco),     INTENT(IN)     :: noco
      TYPE(t_potden),   INTENT(IN)     :: pot
      TYPE(t_greensf),  INTENT(IN)     :: gdft(:) !green's function calculated from the Kohn-Sham system
      TYPE(t_hub1data), INTENT(INOUT)  :: hub1data
      TYPE(t_results),  INTENT(INOUT)  :: results
      TYPE(t_potden),   INTENT(INOUT)  :: den

      INTEGER :: i_hia,nType,l,occDFT_INT,ispin,m,i_exc,n
      INTEGER :: io_error,ierr
      INTEGER :: indStart,indEnd
      INTEGER :: hubbardioUnit
      INTEGER :: n_hia_task,extra,i_hia_start,i_hia_end
      REAL    :: mu_dc
      LOGICAL :: l_firstIT_HIA,l_ccfexist,l_bathexist

      CHARACTER(len=300) :: cwd,path,folder,xPath
      TYPE(t_greensf),ALLOCATABLE :: gu(:)
      TYPE(t_selfen), ALLOCATABLE :: selfen(:)

#ifdef CPP_HDF
      INTEGER(HID_T)     :: greensf_fileID
#endif

      REAL    :: f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL    :: f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      REAL    :: occDFT(atoms%n_hia,input%jspins)
      COMPLEX :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,3)
      COMPLEX, ALLOCATABLE :: e(:)
      COMPLEX, ALLOCATABLE :: ctmp(:)

      !Check if the EDsolver library is linked
#ifndef CPP_EDSOLVER
      CALL juDFT_error("No solver linked for Hubbard 1", hint="Link the edsolver library",calledby="hubbard1_setup")
#endif

      IF(mpi%irank.EQ.0) THEN
         !-------------------------------------------
         ! Create the Input for the Hubbard 1 Solver
         !-------------------------------------------

         !Get the working directory
         CALL get_environment_variable('PWD',cwd)
         path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(hubbard1CalcFolder))
         CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)))
         !Remove everything from the last iteration (Good Idea??)
         CALL SYSTEM('rm -rf ' // TRIM(ADJUSTL(path)) // "/*")

         !Positions of the DFT+HIA elements in all DFT+U related arrays
         indStart = atoms%n_u+1
         indEnd   = atoms%n_u+atoms%n_hia

         ! calculate slater integrals from u and j
         CALL uj2f(input%jspins,atoms%lda_u(indStart:indEnd),atoms%n_hia,f0,f2,f4,f6)

         DO ispin = 1, 1 ! input%jspins
            f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
            f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
            f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
            f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
         END DO

         DO i_hia = 1, atoms%n_hia

            l = atoms%lda_u(atoms%n_u+i_hia)%l
            nType = atoms%lda_u(atoms%n_u+i_hia)%atomType

            IF(ALL(ABS(gdft(i_hia)%gmmpMat).LT.1e-12)) THEN
               CALL juDFT_error("Hubbard-1 has no DFT greensf available",calledby="hubbard1_setup")
            ENDIF

            !Create Subfolder (if there are multiple Hubbard 1 procedures)
            CALL hubbard1_path(atoms,i_hia,folder)
            WRITE(xPath,*) TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(xPath)))

            !-------------------------------------------------------
            ! Calculate the DFT occupation of the correlated shell
            !-------------------------------------------------------
            CALL occmtx(gdft(i_hia),gfinp,input,mmpMat(:,:,i_hia,:))

            !For the first iteration we can fix the occupation and magnetic moments in the inp.xml file
            l_firstIT_HIA = hub1data%iter.EQ.1 .AND.ALL(ABS(den%mmpMat(:,:,indStart:indEnd,:)).LT.1e-12)
            IF(l_firstIT_HIA) THEN
               IF(hub1inp%init_occ(i_hia) > -9e98) THEN
                  occDFT(i_hia,:) = hub1inp%init_occ(i_hia)/input%jspins
               ELSE
                  occDFT(i_hia,:) = 0.0
                  DO ispin = 1, input%jspins
                     DO m = -l, l
                        occDFT(i_hia,ispin) = occDFT(i_hia,ispin) + REAL(mmpMat(m,m,i_hia,ispin))
                     ENDDO
                  ENDDO
               ENDIF

               DO i_exc = 1, hub1inp%n_exc(i_hia)
                  IF(hub1inp%init_mom(i_hia,i_exc) > -9e98) THEN
                     hub1data%mag_mom(i_hia,i_exc) = hub1inp%init_mom(i_hia,i_exc)
                  ENDIF
               ENDDO
            ELSE
               occDFT(i_hia,:) = 0.0
               DO ispin = 1, input%jspins
                  DO m = -l, l
                     occDFT(i_hia,ispin) = occDFT(i_hia,ispin) + REAL(mmpMat(m,m,i_hia,ispin))
                  ENDDO
               ENDDO
            ENDIF
            !Nearest Integer occupation
            occDFT_INT = ANINT(SUM(occDFT(i_hia,:)))

            !Initial Information (We are already on irank 0)
            WRITE(oUnit,*)
            WRITE(oUnit,9010) nType
9010        FORMAT("Setup for Hubbard 1 solver for atom ", I3, ": ")
            WRITE(oUnit,"(A)") "Everything related to the solver (e.g. mu_dc) is given in eV"
            WRITE(oUnit,*)
            WRITE(oUnit,"(A)") "Occupation from DFT-Green's function:"
            WRITE(oUnit,9020) 'spin-up','spin-dn'
9020        FORMAT(TR8,A7,TR3,A7)
            WRITE(oUnit,9030) occDFT(i_hia,:)
9030        FORMAT(TR7,f8.4,TR2,f8.4)

            !--------------------------------------------------------------------------
            ! Calculate the chemical potential for the solver
            ! This is equal to the double-counting correction used in DFT+U
            !--------------------------------------------------------------------------
            ! V_FLL = U (n - 1/2) - J (n - 1) / 2
            ! V_AMF = U n/2 + 2l/[2(2l+1)] (U-J) n
            !--------------------------------------------------------------------------
            CALL mudc(atoms%lda_u(atoms%n_u+i_hia),occDFT(i_hia,:),input%jspins,mu_dc)

            !-------------------------------------------------------
            ! Check for additional input files
            !-------------------------------------------------------
            !Is a crystal field matrix present in the work directory (overwrites the calculated matrix)
            INQUIRE(file=TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(cfg_file_ccf)),exist=l_ccfexist)
            IF(l_ccfexist) CALL read_ccfmat(TRIM(ADJUSTL(cwd)),hub1data%ccfmat(i_hia,-l:l,-l:l),l)
            !Is a bath parameter file present
            INQUIRE(file=TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(cfg_file_bath)),exist=l_bathexist)
            !Copy the bath file to the Hubbard 1 solver if its present
            IF(l_bathexist) CALL SYSTEM('cp ' // TRIM(ADJUSTL(cfg_file_bath)) // ' ' // TRIM(ADJUSTL(xPath)))

            !-------------------------------------------------------
            ! Write the main config files
            !-------------------------------------------------------
            CALL write_hubbard1_input(xPath,i_hia,l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),&
                                      hub1inp,hub1data,mu_dc,occDFT_INT,l_bathexist,l_firstIT_HIA)
         ENDDO
      ENDIF !mpi%irank == 0

      IF(mpi%irank.EQ.0) THEN
         WRITE(*,*) "Calculating new density matrix ..."
      ENDIF

      ALLOCATE(gu(atoms%n_hia))
      ALLOCATE(selfen(atoms%n_hia))

      !Argument order different because occDFT is not allocatable
      CALL mpi_bc(mpi%irank,mpi%mpi_comm,occDFT)
      !Broadcast important stuff
      DO i_hia = 1, atoms%n_hia
         CALL gu(i_hia)%init(gdft(i_hia)%elem,gfinp,input,contour_in=gdft(i_hia)%contour)
         CALL selfen(i_hia)%init(lmaxU_const,gdft(i_hia)%contour%nz)
      ENDDO

#ifdef CPP_MPI
      !distribute the individual hubbard1 elements over the ranks
      n_hia_task = FLOOR(REAL(atoms%n_hia)/(mpi%isize))
      extra = atoms%n_hia - n_hia_task*mpi%isize
      i_hia_start = mpi%irank*n_hia_task + 1 + extra
      i_hia_end   =(mpi%irank+1)*n_hia_task   + extra
      IF(mpi%irank < extra) THEN
         i_hia_start = i_hia_start - (extra - mpi%irank)
         i_hia_end   = i_hia_end   - (extra - mpi%irank - 1)
      ENDIF
#else
      i_hia_start = 1
      i_hia_end   = atoms%n_hia
#endif

#ifdef CPP_MPI
      !Make sure that the ranks are synchronized
      CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

      mmpMat = cmplx_0
      !------------------------------------------------------------
      ! This loop runs the solver
      !------------------------------------------------------------
      DO i_hia = i_hia_start, i_hia_end

         IF(i_hia > atoms%n_hia .OR. i_hia < 1) CYCLE

         nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
         l = atoms%lda_u(atoms%n_u+i_hia)%l

         CALL get_environment_variable('PWD',cwd)
         CALL hubbard1_path(atoms,i_hia,folder)
         WRITE(xPath,*) TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))

         ALLOCATE(e(gdft(i_hia)%contour%nz),source=cmplx_0)

         CALL timestart("Hubbard 1: EDsolver")
         !We have to change into the Hubbard1 directory so that the solver routines can read the config
         CALL CHDIR(TRIM(ADJUSTL(xPath)))
#ifdef CPP_EDSOLVER
         !Open the output file for the solver
         hubbardioUnit = 4000+i_hia
         OPEN(unit=hubbardioUnit, file=TRIM(ADJUSTL(xPath)) // TRIM(ADJUSTL(hubbard1Outfile)),&
              status="replace", action="write", iostat=io_error)
         IF(io_error/=0) CALL juDFT_error("Error in opening EDsolver out file",calledby="hubbard1_setup")
         e = gdft(i_hia)%contour%e*hartree_to_ev_const
         CALL EDsolver_from_cfg(2*(2*l+1),gdft(i_hia)%contour%nz,e,selfen(i_hia)%data(:,:,:,1),1,hubbardioUnit)
         !---------------------------------------------------
         ! Calculate selfenergy on lower contour explicitly
         ! Mainly out of paranoia :D
         ! No rediagonalization (last argument switches this)
         !---------------------------------------------------
         e = conjg(gdft(i_hia)%contour%e)*hartree_to_ev_const
         CALL EDsolver_from_cfg(2*(2*l+1),gdft(i_hia)%contour%nz,e,selfen(i_hia)%data(:,:,:,2),0,hubbardioUnit)
         CLOSE(hubbardioUnit, iostat=io_error)
         IF(io_error/=0) CALL juDFT_error("Error in closing EDsolver out file",calledby="hubbard1_setup")
#endif
         CALL CHDIR(TRIM(ADJUSTL(cwd)))
         CALL timestop("Hubbard 1: EDsolver")

         DEALLOCATE(e)

         !-------------------------------------------
         ! Postprocess selfenergy
         !-------------------------------------------
         CALL selfen(i_hia)%postProcess(input%jspins,noco%l_mtNocoPot.AND.gfinp%l_mperp,pot%mmpMat(:,:,atoms%n_u+i_hia,:))

         !----------------------------------------------------------------------
         ! Solution of the Dyson Equation
         !----------------------------------------------------------------------
         ! G(z)^(-1) = G_0(z)^(-1) - mu - Sigma(z)
         !----------------------------------------------------------------------
         ! Sigma(z) is the self-energy from the impurity solver
         ! We introduce an additional chemical potential mu, which is determined
         ! so that the occupation of the correlated orbital does not change
         !----------------------------------------------------------------------
         CALL timestart("Hubbard 1: Add Selfenergy")
         CALL add_selfen(gdft(i_hia),i_hia,selfen(i_hia),gfinp,input,noco,&
                         occDFT(i_hia,:),gu(i_hia),mmpMat(:,:,i_hia,:))
         CALL timestop("Hubbard 1: Add Selfenergy")

      ENDDO

      !Collect the impurity Green's Function
      DO i_hia = 1, atoms%n_hia
         CALL gu(i_hia)%collect(mpi%mpi_comm)
         CALL selfen(i_hia)%collect(mpi%mpi_comm)
      ENDDO


#ifdef CPP_HDF
      IF(mpi%irank.EQ.0) THEN
         !------------------------------
         !Write out DFT Green's Function
         !------------------------------
         CALL timestart("Hubbard 1: IO/Write")
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms, inFilename="greensf_DFT.hdf")
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, &
                               GREENSF_HUBBARD_CONST, gdft, mmpmat)
         CALL closeGreensFFile(greensf_fileID)

         !-------------------------------------
         !Write out correlated Green's Function
         !-------------------------------------
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms, inFilename="greensf_IMP.hdf")
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, &
                              GREENSF_HUBBARD_CONST, gu, mmpmat,selfen=selfen)
         CALL closeGreensFFile(greensf_fileID)
         CALL timestop("Hubbard 1: IO/Write")
      ENDIF
#endif


#ifdef CPP_MPI
      !Collect the density matrix to rank 0
      n = SIZE(mmpMat)
      ALLOCATE(ctmp(n))
      CALL MPI_REDUCE(mmpMat,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF(mpi%irank.EQ.0) CALL CPP_BLAS_ccopy(n,ctmp,1,mmpMat,1)
      DEALLOCATE(ctmp)
#endif

      !--------------------------------------------------------------------
      ! Calculate Distances from last density matrix and update den%mmpmat
      !--------------------------------------------------------------------
      results%last_mmpMatdistance = 0.0
      results%last_occdistance = 0.0

      IF(mpi%irank.EQ.0) THEN
         DO i_hia = 1, atoms%n_hia
            CALL hubbard1Distance(den%mmpMat(:,:,atoms%n_u+i_hia,:),mmpMat(:,:,i_hia,:),input,gfinp,results)
            DO ispin = 1, MERGE(3,input%jspins,gfinp%l_mperp)
               den%mmpMat(-lmaxU_const:,-lmaxU_const:,atoms%n_u+i_hia,ispin) = mmpMat(-lmaxU_const:,-lmaxU_const:,i_hia,ispin)
            ENDDO
         ENDDO
      ENDIF

      !Broadcast the density matrix
      CALL mpi_bc(den%mmpMat,mpi%irank,mpi%mpi_comm)

   END SUBROUTINE hubbard1_setup

   SUBROUTINE hubbard1_path(atoms,i_hia,xPath)

      !Defines the folder structure
      ! The Solver is run in the subdirectories
      ! Hubbard1/ if only one Hubbard1 procedure is run
      ! Hubbard1/atom_label_l if there are more

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      INTEGER,             INTENT(IN)  :: i_hia
      CHARACTER(len=300),  INTENT(OUT) :: xPath

      CHARACTER(len=300) :: folder,fmt
      CHARACTER(len=1),PARAMETER :: spdfg(0:4) = ['s','p','d','f','g']
      INTEGER nType,l

      nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
      l = atoms%lda_u(atoms%n_u+i_hia)%l
      xPath = TRIM(ADJUSTL(hubbard1CalcFolder))
      IF(atoms%n_hia>1) THEN
         WRITE(fmt,'("(A",I2.2,",A1,A1,A1)")') LEN(TRIM(ADJUSTL(atoms%label(nType))))
         WRITE(folder,fmt) TRIM(ADJUSTL(atoms%label(nType))),"_",spdfg(l),"/"
      ELSE
         folder=""
      ENDIF
      xPath = TRIM(ADJUSTL(xPath)) // "/" // folder

   END SUBROUTINE hubbard1_path

END MODULE m_hubbard1_setup
