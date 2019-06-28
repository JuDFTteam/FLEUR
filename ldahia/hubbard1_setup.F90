MODULE m_hubbard1_setup

   USE m_juDFT

   IMPLICIT NONE

   LOGICAL, PARAMETER :: l_debug = .TRUE.  !Enable/Disable Debug outputs like dependency of occupation on chemical potential shift 

   CONTAINS

   SUBROUTINE hubbard1_setup(iterHIA,atoms,hub1,sym,mpi,noco,input,usdus,pot,gdft,results)

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
#ifdef CPP_MPI
      INCLUDE "mpif.h"
#endif

      INTEGER,          INTENT(INOUT)  :: iterHIA !number of iteration 
      TYPE(t_atoms),    INTENT(IN)     :: atoms
      TYPE(t_hub1ham),  INTENT(IN)     :: hub1
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

      INTEGER i_hia,i_gf,n,l,n_occ,ispin,m,iz,k,j
      INTEGER io_error,ierr
      REAL    mu_dc,e_lda_hia

      CHARACTER(len=300) :: cwd,path,folder,message
      CHARACTER(len=8)   :: l_type*2,l_form*9

      TYPE(t_greensf)    :: gu 

      REAL     f0(atoms%n_hia,input%jspins),f2(atoms%n_hia,input%jspins)
      REAL     f4(atoms%n_hia,input%jspins),f6(atoms%n_hia,input%jspins)
      REAL     u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
               -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia)

      COMPLEX  mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX  mmpMat_in(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      COMPLEX  selfen(atoms%n_hia,gdft%nz,2*(2*lmaxU_const+1),2*(2*lmaxU_const+1))
      REAL     n_l(atoms%n_hia,input%jspins)
      LOGICAL  l_exist

      !TODO: We don't need to calculate the green's function in every iteration

      !Read in density matrix from n_mmpmat_hubbard1
      !We use a separate file because we need to keep track of the distances and the iteration between different runs of fleur
      INQUIRE(file="n_mmpmat_hubbard1",exist=l_exist)
      IF(l_exist) THEN
         OPEN(unit = 1337, file="n_mmpmat_hubbard1",status="old",form="formatted",action="read")
         READ(1337,9110) iterHIA,results%last_occdistance,results%last_mmpMatdistance
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

         CALL v_mmp(sym,atoms,atoms%lda_hia(:),atoms%n_hia,input%jspins,.true.,mmpMat,&
         u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),e_lda_hia)

         IF(hub1%l_runthisiter.AND.ANY(gdft%gmmpMat(:,:,:,:,:,:).NE.0.0)) THEN 
            !IF(gdft%mode.EQ.2) CALL juDFT_error("This energy contour is not supported at the moment for DFT+Hubbard1",calledby="hubbard1_setup")
            !The onsite green's function was calculated but the solver 
            !was not yet run
            !--> write out the configuration for the hubbard 1 solver 
            CALL gu%init(input,lmaxU_const,atoms,.true.,noco,nz_in=gdft%nz, e_in=gdft%e,de_in=gdft%de,matsub_in=gdft%nmatsub)

            !Get the working directory
            CALL get_environment_variable('PWD',cwd)
            !Create a folder where to store the files from the solver
            !Is this applicable on all systems where fleur can be run?
            iterHIA = iterHIA + 1
            WRITE(folder,"(A5,I3.3)") "hub1_" , iterHIA
            path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
            CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)))

            DO i_hia = 1, atoms%n_hia
               n = atoms%lda_hia(i_hia)%atomType
               l = atoms%lda_hia(i_hia)%l
               CALL indexgf(atoms,l,n,i_gf)
               !Create a subdirectory for the atomType and shell
               WRITE(folder,"(A4,I2.2,A2,I1.1)") "atom",n,"_l",l
               CALL SYSTEM('mkdir -p ' // TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)))
            
               !calculate the occupation of the correlated shell
               CALL occmtx(gdft,i_gf,atoms,sym,input,results%ef,mmpMat(:,:,i_hia,:))
               n_l(i_hia,:) = 0.0
               DO ispin = 1, input%jspins
                  DO m = -l, l
                     n_l(i_hia,ispin) = n_l(i_hia,ispin) + mmpMat(m,m,i_hia,ispin)
                  ENDDO
               ENDDO
               WRITE(6,*)
               WRITE(6,9010) n
               WRITE(6,"(A)") "Everything related to the solver (e.g. mu_dc) is given in eV"
               WRITE(6,*)
               WRITE(6,"(A)") "Occupation from DFT-Green's function:"
               WRITE(6,9020) 'spin-up','spin-dn'
               WRITE(6,9030) n_l(i_hia,:)
               !calculate the chemical potential
               CALL mudc(atoms%lda_hia(i_hia)%U,atoms%lda_hia(i_hia)%J,l,n_l(i_hia,:),atoms%lda_hia(i_hia)%l_amf,mu_dc,input%jspins)
               
               n_occ = ANINT(SUM(n_l(i_hia,:)))
               n_l(i_hia,1) = SUM(n_l(i_hia,:))
               !Check wether the hubbard 1 solver was run:
               INQUIRE(file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/" // "se.atom",exist=l_exist)
               IF(mpi%irank.EQ.0) THEN
                  IF(l_exist) THEN
                     CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,1:gdft%nz-gdft%nmatsub,1:2*(2*l+1),1:2*(2*l+1)),gdft%nz-gdft%nmatsub,2*(2*l+1),.false.)
                     IF(gdft%nmatsub.GT.0) CALL read_selfen(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",selfen(i_hia,gdft%nz-gdft%nmatsub:gdft%nz,1:2*(2*l+1),1:2*(2*l+1)),gdft%nmatsub,2*(2*l+1),.true.)
                  ELSE
                     CALL write_hubbard1_input(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",l,f0(i_hia,1),f2(i_hia,1),f4(i_hia,1),f6(i_hia,1),&
                                                hub1%xi(i_hia),-hub1%bz(i_hia),MAX(1,n_occ-hub1%n_exc),MIN(2*(2*l+1),n_occ+hub1%n_exc),hub1%beta,mu_dc,hub1%l_ccf(i_hia),&
                                                gdft%nz-gdft%nmatsub,gdft%nmatsub,REAL(gdft%e(1)),REAL(gdft%e(gdft%nz-gdft%nmatsub)),AIMAG(gdft%e(1)))
                     
                     !WRITE(*,"(7f14.8)") hub1%ccfmat(i_hia,:,:)
                     IF(hub1%l_ccf(i_hia)) THEN
                        CALL write_ccfmat(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/",hub1%ccfmat(i_hia,-l:l,-l:l),l)
                     ENDIF

                     IF(gdft%mode.NE.1) THEN
                        OPEN(unit=1337,file=TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/" // "contour.dat",status="replace",action="write")
                        DO iz = 1, gdft%nz
                           WRITE(1337,"(2f14.8)") REAL(gdft%e(iz))*hartree_to_ev_const, AIMAG(gdft%e(iz))*hartree_to_ev_const
                        ENDDO
                        CLOSE(unit= 1337)
                     ENDIF     
                  ENDIF
               ENDIF
            ENDDO
            IF(mpi%irank.EQ.0) CALL CHDIR(TRIM(ADJUSTL(cwd)))

            IF(l_exist) THEN
               CALL add_selfen(gdft,gu,selfen,atoms,hub1,sym,input,results%ef,n_l(:,1),mu_dc/hartree_to_ev_const,&
                              pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),mmpMat)
               ! calculate potential matrix and total energy correction
               CALL v_mmp(sym,atoms,atoms%lda_hia,atoms%n_hia,input%jspins,.true.,mmpMat,&
                     u,f0,f2,pot%mmpMat(:,:,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),results%e_ldau)
               !
               ! Output the density of states from the two green's functions (DEBUG)
               !
               IF(l_debug) THEN
                  DO i_hia = 1, atoms%n_hia
                     n = atoms%lda_hia(i_hia)%atomType
                     l = atoms%lda_hia(i_hia)%l
                     WRITE(folder,"(A5,I3.3)") "hub1_" , iterHIA
                     path = TRIM(ADJUSTL(cwd)) // "/" // TRIM(ADJUSTL(folder))
                     WRITE(folder,"(A4,I2.2,A2,I1.1)") "atom",n,"_l",l
                     CALL CHDIR(TRIM(ADJUSTL(path)) // "/" // TRIM(ADJUSTL(folder)) // "/")
                     CALL indexgf(atoms,l,n,i_gf)
                     CALL ldosmtx("gdft",gdft,i_gf,atoms,sym,input)
                     CALL ldosmtx("g",gu,i_gf,atoms,sym,input)
                  ENDDO
               ENDIF
               !
               !Read in the last density matrix
               !
               CALL CHDIR(TRIM(ADJUSTL(cwd)))
               INQUIRE(file="n_mmpmat_hubbard1",exist = l_exist)
               IF(l_exist) THEN
                  OPEN(unit=1337,file="n_mmpmat_hubbard1",status="old",action="read",iostat=io_error)
                  IF(io_error.NE.0) CALL juDFT_error("IO-error in density matrix",calledby="hubbard1_setup")
                  READ(1337,9110)
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
                     WRITE(1337,9110) iterHIA,results%last_occdistance,results%last_mmpMatdistance
                  WRITE(1337,"(7f14.8)") mmpMat(:,:,:,:)
                  CLOSE(unit=1337)
                  WRITE(*,*) "Hubbard 1 Iteration: ", iterHIA
                  WRITE(*,*)  "Occ. Distance: ", results%last_occdistance, &
                              "Mat. Distance: ", results%last_mmpMatdistance
               ENDIF
            ELSE
               WRITE(message,9100) iterHIA
               CALL juDFT_END(TRIM(ADJUSTL(message)),mpi%irank)
            ENDIF
         ENDIF
#ifdef CPP_MPI
         CALL MPI_BCAST(pot%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u+1:atoms%n_u+atoms%n_hia,:),&
                        49*atoms%n_hia*input%jspins,MPI_DOUBLE_COMPLEX,0,mpi%mpi_comm,ierr)
#endif 
         !Write out the density matrix and potential matrix (compare u_setup.f90)
         IF (mpi%irank.EQ.0) THEN
            DO ispin = 1,input%jspins
               WRITE (6,'(a7,i3)') 'spin #',ispin
               DO i_hia = 1, atoms%n_hia
                  n = atoms%lda_hia(i_hia)%atomType
                  l = atoms%lda_hia(i_hia)%l
                  WRITE (l_type,'(i2)') 2*(2*l+1)
                  l_form = '('//l_type//'f12.7)'
                  WRITE (6,'(a20,i3)') 'n-matrix for atom # ',n
                  WRITE (6,l_form) ((mmpMat(k,j,atoms%n_u+i_hia,ispin),k=-l,l),j=-l,l)
                  WRITE (6,'(a20,i3)') 'V-matrix for atom # ',n
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

      ELSE
         !occupation matrix is zero and LDA+Hubbard 1 shouldn't be run yet
         !There is nothing to be done yet just set the potential correction to 0
         pot%mmpMat(:,:,atoms%n_u+1:atoms%n_hia+atoms%n_u,:) = CMPLX(0.0,0.0)
      ENDIF
      
      !FORMAT Statements:
9010  FORMAT("Setup for Hubbard 1 solver for atom ", I3, ": ")
9020  FORMAT(TR8,A7,TR3,A7)
9030  FORMAT(TR7,f8.4,TR2,f8.4)
9100  FORMAT("Hubbard 1 input written into hub1_",I3.3)
9110  FORMAT("Hubbard 1 Iteration ", I3, " Last Distance in Occupation: ",f14.8," Last density matrix Distance: ", f14.8)

   END SUBROUTINE hubbard1_setup

END MODULE m_hubbard1_setup