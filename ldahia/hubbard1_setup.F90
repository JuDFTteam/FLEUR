MODULE m_hubbard1_setup

   USE m_juDFT

   IMPLICIT NONE

   LOGICAL, PARAMETER :: l_debug = .TRUE.  !Enable/Disable Debug outputs like dependency of occupation on chemical potential shift 
   CHARACTER(len=30), PARAMETER :: main_folder = "Hubbard1"



   CONTAINS

   SUBROUTINE hubbard1_setup(atoms,hub1,sym,mpi,noco,input,usdus,pot,gdft,results)

      USE m_types
      USE m_hubbard1_io
      USE m_uj2f
      USE m_constants
      USE m_gfcalc
      USE m_umtx
      USE m_vmmp
      USE m_add_selfen
      USE m_mudc
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
      TYPE(t_potden),   INTENT(INOUT)  :: pot
      TYPE(t_results),  INTENT(INOUT)  :: results
      TYPE(t_greensf),  INTENT(IN)     :: gdft !green's function calculated from the Kohn-Sham system
      
#ifdef CPP_MPI
      EXTERNAL MPI_BCAST
#endif

      INTEGER i_hia,nType,l,n_occ,ispin,m,iz,k,j,i_exc
      INTEGER io_error,ierr
      REAL    mu_dc,e_lda_hia,exc

      CHARACTER(len=300) :: cwd,path,folder,message
      CHARACTER(len=8)   :: l_type*2,l_form*9
      CHARACTER(len=1)   :: l_name(0:3)
      CHARACTER(len=30)  :: fmt

      TYPE(t_greensf)    :: gu 

      REAL     f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL     f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      REAL     u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
               -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia)

      COMPLEX  mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX  mmpMat_in(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX  selfen(atoms%n_hia,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1),gdft%nz)
      REAL     n_l(atoms%n_hia,input%jspins)
      LOGICAL  l_selfenexist,l_exist,l_linkedsolver

      l_name(0:3) = (/"s","p","d","f"/)

      !TODO: We don't need to calculate the green's function in every iteration

      !Read in density matrix from n_mmpmat_hubbard1
      !We use a separate file because we need to keep track of the distances and the iteration between different runs of fleur
      INQUIRE(file="n_mmpmat_hubbard1",exist=l_exist)
      IF(l_exist) THEN
         OPEN(unit = 1337, file="n_mmpmat_hubbard1",status="old",form="formatted",action="read")
         READ(1337,"(7f14.8)") mmpMat(:,:,:,:)
         CLOSE(unit=1337)
      ELSE
        mmpMat(:,:,:,:) = 0.0
      ENDIF

      IF(ANY(mmpMat(:,:,:,:).NE.0.0).OR.hub1%l_runthisiter) THEN
         !Get the slater integrals from the U and J parameters
         CALL uj2f(input%jspins,atoms%lda_hia(:),atoms%n_hia,f0,f2,f4,f6)
         f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
         f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
         f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
         f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
         CALL umtx(atoms%lda_hia(:),atoms%n_hia,f0(:,1),f2(:,1),f4(:,1),f6(:,1),u)

         CALL v_mmp(sym,atoms,atoms%lda_hia(:),atoms%n_hia,input%jspins,.FALSE.,mmpMat,&
         u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),e_lda_hia)

         IF(hub1%l_runthisiter.AND.(ANY(gdft%gmmpMat(:,:,:,:,:,:).NE.0.0))) THEN 
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
               nType = atoms%lda_hia(i_hia)%atomType
               l = atoms%lda_hia(i_hia)%l

               !Create a subdirectory for the atomType and shell
               IF(atoms%n_hia>1) THEN
                  WRITE(fmt,'("(A",I2.2,",A1,A1,A1")') LEN(TRIM(ADJUSTL(atoms%label(nType))))
                  WRITE(folder,fmt) atoms%label(nType),"_",l_name(l),"/"
                  CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
               ELSE 
                  folder=""
               ENDIF
               IF(hub1%iter.EQ.1.AND.ALL(mmpMat.EQ.0.0)) THEN
                  n_l(i_hia,:) = hub1%init_occ(i_hia)/input%jspins
                  DO i_exc = 1, hub1%n_exc_given(i_hia)
                     hub1%mag_mom(i_hia,i_exc) = hub1%init_mom(i_hia,i_exc)
                  ENDDo
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

               WRITE(6,*)
               WRITE(6,9010) nType
               WRITE(6,"(A)") "Everything related to the solver (e.g. mu_dc) is given in eV"
               WRITE(6,*)
               WRITE(6,"(A)") "Occupation from DFT-Green's function:"
               WRITE(6,9020) 'spin-up','spin-dn'
               WRITE(6,9030) n_l(i_hia,:)

               !calculate the chemical potential
               CALL mudc(atoms%lda_hia(i_hia)%U,atoms%lda_hia(i_hia)%J,l,n_l(i_hia,:),atoms%lda_hia(i_hia)%l_amf,mu_dc,input%jspins)
               n_occ = ANINT(SUM(n_l(i_hia,:)))

               !Check wether the hubbard 1 solver was run:(old version)
               INQUIRE(file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "se.atom",exist=l_selfenexist)
               IF(mpi%irank.EQ.0.AND..NOT.l_selfenexist) THEN
#ifdef CPP_EDSOLVER
                  CALL hubbard1_input(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",i_hia,l,&
                                            f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),hub1,mu_dc,n_occ,.true.)
#else
                  CALL hubbard1_input(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",i_hia,l,&
                                            f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),hub1,mu_dc,n_occ,.false.)
                  !If no Solver is linked we assume that the old solver is used and we write out some additional files
                  IF(hub1%ccf(i_hia).NE.0.0) THEN
                     CALL write_ccfmat(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",hub1%ccfmat(i_hia,-l:l,-l:l),l)
                  ENDIF
                  IF(gdft%mode.NE.1) THEN
                     OPEN(unit=1337,file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/" // "contour.dat",status="replace",action="write")
                     DO iz = 1, gdft%nz
                        WRITE(1337,"(2f14.8)") REAL(gdft%e(iz))*hartree_to_ev_const, AIMAG(gdft%e(iz))*hartree_to_ev_const
                     ENDDO
                     CLOSE(unit= 1337)
                  ENDIF   
#endif
               ENDIF
            ENDDO

            
            DO i_hia = 1, atoms%n_hia 
               nType = atoms%lda_hia(i_hia)%atomType
               l = atoms%lda_hia(i_hia)%l
               IF(atoms%n_hia>1) THEN
                  WRITE(fmt,'("(A",I2.2,",A1,A1,A1")') LEN(TRIM(ADJUSTL(atoms%label(nType))))
                  WRITE(folder,fmt) atoms%label(nType),"_",l_name(l),"/"
                  CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
               ELSE
                  folder=""
               ENDIF
               INQUIRE(file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "se.atom",exist=l_selfenexist)
               IF(.NOT.l_selfenexist.AND.mpi%irank.EQ.0) THEN
#ifdef CPP_EDSOLVER
                  l_linkedsolver=.TRUE.
                  CALL timestart("Hubbard 1: EDsolver")
                  CALL CHDIR(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
                  CALL EDsolver_from_cfg(2*(2*l+1),gdft%nz,gdft%e(1:gdft%nz)*hartree_to_ev_const,selfen(i_hia,:,:,1:gdft%nz),1)
                  selfen(i_hia,:,:,:) = selfen(i_hia,:,:,:)/hartree_to_ev_const
                  CALL CHDIR(cwd)
                  CALL timestop("Hubbard 1: EDsolver")
#else 
                  CALL juDFT_END("Hubbard1 input has been written into Hubbard1/ (No Solver linked)",mpi%irank)
#endif
               ENDIF
            ENDDO

            IF(l_selfenexist.OR.l_linkedsolver) THEN
               CALL timestart("Hubbard 1: Add Selfenenergy")
               IF(l_selfenexist) CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",&
                                                selfen(i_hia,1:2*(2*l+1),1:2*(2*l+1),1:gdft%nz),gdft%nz,2*(2*l+1),.false.)
               CALL add_selfen(gdft,gu,selfen,atoms,noco,hub1,sym,input,results%ef,n_l,mu_dc/hartree_to_ev_const,&
                              pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),mmpMat)
               CALL timestop("Hubbard 1: Add Selfenenergy")
               ! calculate potential matrix and total energy correction
               CALL v_mmp(sym,atoms,atoms%lda_hia,atoms%n_hia,input%jspins,.FALSE.,mmpMat,&
                     u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),results%e_ldau)
               !
               !Read in the last density matrix
               !
               INQUIRE(file="n_mmpmat_hubbard1",exist = l_exist)
               IF(l_exist) THEN
                  OPEN(unit=1337,file="n_mmpmat_hubbard1",status="old",action="read",iostat=io_error)
                  IF(io_error.NE.0) CALL juDFT_error("IO-error in density matrix",calledby="hubbard1_setup")
                  READ(1337,"(7f14.8)") mmpMat_in(:,:,:,:)
                  CLOSE(unit=1337)
               ELSE
                  mmpMat_in = 0.0
               ENDIF
               CALL n_mmp_dist(mmpMat_in,mmpMat,atoms%n_hia,results,input%jspins)
               IF(mpi%irank.EQ.0) THEN
                  !Write out the density matrix and the additional inforamtion (current iteration, distances)
                  OPEN(unit=1337,file="n_mmpmat_hubbard1",status="replace",action="write",iostat=io_error)
                  IF(io_error.NE.0) CALL juDFT_error("IO-error in density matrix",calledby="hubbard1_setup")
                  WRITE(1337,"(7f14.8)") mmpMat(:,:,:,:)
                  CLOSE(unit=1337)
               ENDIF
            ENDIF
            !Write out the density matrix and potential matrix (compare u_setup.f90)
            IF (mpi%irank.EQ.0) THEN
               DO ispin = 1,input%jspins
                  WRITE (6,'(a7,i3)') 'spin #',ispin
                  DO i_hia = 1, atoms%n_hia
                     nType = atoms%lda_hia(i_hia)%atomType
                     l = atoms%lda_hia(i_hia)%l
                     WRITE (l_type,'(i2)') 2*(2*l+1)
                     l_form = '('//l_type//'f12.7)'
                     WRITE (6,'(a20,i3)') 'n-matrix for atom # ',nType
                     WRITE (6,l_form) ((mmpMat(k,j,atoms%n_u+i_hia,ispin),k=-l,l),j=-l,l)
                     WRITE (6,'(a20,i3)') 'V-matrix for atom # ',nType
                     IF (atoms%lda_hia(i_hia)%l_amf) THEN
                        WRITE (6,*) 'using the around-mean-field limit '
                     ELSE
                        WRITE (6,*) 'using the atomic limit of LDA+U '
                     ENDIF
                     WRITE (6,l_form) ((pot%mmpMat(k,j,atoms%n_u+i_hia,ispin),k=-l,l),j=-l,l)
                  END DO
               END DO
               WRITE (6,*) results%e_ldau
            ENDIF
         ENDIF
      ELSE IF(mpi%irank.NE.0) THEN
         pot%mmpMat(:,:,atoms%n_u+1:atoms%n_hia+atoms%n_u,:) = CMPLX(0.0,0.0)
         !If we are on a different mpi%irank and no solver is linked we need to call juDFT_end here if the solver was not run
         !kind of a weird workaround (replace with something better)
#ifdef CPP_EDSOLVER 
         !Do nothing and go to the MPI_BCAST
#else 
         CALL get_environment_variable('PWD',cwd)
         path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(main_folder))
         DO i_hia = 1, atoms%n_hia 
            nType = atoms%lda_hia(i_hia)%atomType
            l = atoms%lda_hia(i_hia)%l
            IF(atoms%n_hia>1) THEN
               WRITE(fmt,'("(A",I2.2,",A1,A1,A1")') LEN(TRIM(ADJUSTL(atoms%label(nType))))
               WRITE(folder,fmt) atoms%label(nType),"_",l_name(l),"/"
               CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
            ELSE
               folder=""
            ENDIF
            INQUIRE(file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "se.atom",exist=l_selfenexist)
            IF(.NOT.l_selfenexist) CALL juDFT_END("Hubbard1 input has been written into Hubbard1/ (No Solver linked)",mpi%irank)
         ENDDO
         !If we are here the solver was run and we go to MPI_BCAST
#endif
      ELSE
         !occupation matrix is zero and LDA+Hubbard 1 shouldn't be run yet
         !There is nothing to be done yet just set the potential correction to 0
         pot%mmpMat(:,:,atoms%n_u+1:atoms%n_hia+atoms%n_u,:) = CMPLX(0.0,0.0)
      ENDIF
#ifdef CPP_MPI
      CALL MPI_BCAST(pot%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),&
                     49*atoms%n_hia*input%jspins,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
#endif 
      
      !FORMAT Statements:
9010  FORMAT("Setup for Hubbard 1 solver for atom ", I3, ": ")
9020  FORMAT(TR8,A7,TR3,A7)
9030  FORMAT(TR7,f8.4,TR2,f8.4)

   END SUBROUTINE hubbard1_setup

END MODULE m_hubbard1_setup