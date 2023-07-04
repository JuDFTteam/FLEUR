!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_init
#ifdef CPP_MPI
   use mpi
#endif
   IMPLICIT NONE
CONTAINS
   SUBROUTINE fleur_init(fmpi, fi, sphhar, stars, nococonv, forcetheo, enpara, xcpot, results, wann, hybdat, mpdata, filename_add)
      USE m_types
      USE m_test_performance
      use m_store_load_hybrid
      USE m_fleurinput_read_xml
      USE m_fleurinput_mpi_bc
      USE m_types_mpinp
      USE m_judft
      USE m_juDFT_init
      USE m_init_wannier_defaults
      USE m_dwigner
      USE m_ylm
      !USE m_InitParallelProcesses
      USE m_xmlOutput
      USE m_constants
      USE m_writeOutParameters
      USE m_setupMPI
      USE m_cdn_io
      USE m_fleur_info
      USE m_mixing_history
      USE m_checks
      USE m_writeOutHeader
      !USE m_fleur_init_old
      USE m_types_xcpot_inbuild
      USE m_make_stars
      USE m_make_sphhar
      USE m_convn
      USE m_efield
      USE m_fleurinput_postprocess
      USE m_make_forcetheo
      USE m_lapwdim
      use m_make_xcpot
      USE m_gaunt, ONLY: gaunt_init
#ifdef CPP_HDF
      USE m_hdf_tools
#endif
      IMPLICIT NONE
      !     Types, these variables contain a lot of data!

      TYPE(t_mpi), INTENT(INOUT):: fmpi
      type(t_fleurinput), intent(out) :: fi
      TYPE(t_sphhar), INTENT(OUT):: sphhar
      TYPE(t_stars), INTENT(OUT):: stars
      TYPE(t_enpara), INTENT(OUT):: enpara
      CLASS(t_xcpot), ALLOCATABLE, INTENT(OUT):: xcpot
      TYPE(t_results), INTENT(OUT):: results
      TYPE(t_wann), INTENT(OUT):: wann
      CLASS(t_forcetheo), ALLOCATABLE, INTENT(OUT)::forcetheo
      TYPE(t_nococonv), INTENT(OUT) :: nococonv
      type(t_hybdat), intent(out) :: hybdat
      type(t_mpdata), intent(out):: mpdata

      CHARACTER(len=100), OPTIONAL, INTENT(IN) :: filename_add

      TYPE(t_enparaXML)::enparaXML
      TYPE(t_forcetheo_data)::forcetheo_data

      TYPE(t_kpts), ALLOCATABLE     :: kptsArray(:)
      INTEGER, ALLOCATABLE          :: xmlElectronStates(:, :)
      INTEGER, ALLOCATABLE          :: atomTypeSpecies(:)
      INTEGER, ALLOCATABLE          :: speciesRepAtomType(:)
      REAL, ALLOCATABLE             :: xmlCoreOccs(:, :, :)
      LOGICAL, ALLOCATABLE          :: xmlPrintCoreStates(:, :)
      !     .. Local Scalars ..
      INTEGER    :: i, n, l, m1, m2, isym, iisym, numSpecies, pc, iAtom, iType, minneigd, outxmlFileID
      INTEGER    :: nbasfcn
      COMPLEX    :: cdum
      CHARACTER(len=4)              :: namex
      CHARACTER(len=12)             :: relcor, tempNumberString
      CHARACTER(LEN=20)             :: filename, tempFilename
      CHARACTER(len=100)            :: filename_add_loc
      CHARACTER(LEN=40)             :: kptsSelection(3)
      CHARACTER(LEN=300)            :: line
      REAL                          :: a1(3), a2(3), a3(3)
      REAL                          :: dtild, phi_add
      LOGICAL                       :: l_found, l_kpts, l_exist, l_krla, l_timeReversalCheck

#ifdef CPP_MPI
      INTEGER ierr(3)
      CALL MPI_COMM_RANK(fmpi%mpi_comm, fmpi%irank, ierr(1))
      CALL MPI_COMM_SIZE(fmpi%mpi_comm, fmpi%isize, ierr(1))
#else
      fmpi%irank = 0; fmpi%isize = 1; fmpi%mpi_comm = 1
#endif
      CALL check_command_line(fmpi)
#ifdef CPP_HDF
      CALL hdf_init()
#endif
      IF (fmpi%irank .EQ. 0) THEN
         filename_add_loc = ""
         IF (PRESENT(filename_add)) filename_add_loc = filename_add
         INQUIRE(file=TRIM(filename_add_loc)//"out.xml", exist=l_exist)
         IF (l_exist) THEN
            tempFilename = "outHistError.xml"
            DO i = 1, 999
               WRITE (tempFilename,'(a,i3.3,a)') 'out-', i, '.xml'
               INQUIRE(file=TRIM(ADJUSTL(tempFilename)), exist=l_found)
               IF (.NOT.l_found) EXIT
            END DO
            IF(.NOT.l_found) THEN
               WRITE(line,'(2a)') 'mv out.xml ', TRIM(ADJUSTL(tempFilename))
               CALL system(TRIM(ADJUSTL(line)))
               WRITE (*,*) 'Moving old out.xml to ', TRIM(ADJUSTL(tempFilename)), '.'
            ELSE
               CALL juDFT_warn("No free out-???.xml file places for storing old out.xml files!")
            END IF
         END IF
         CALL startFleur_XMLOutput(filename_add_loc)
         outxmlFileID = getXMLOutputUnitNumber()
         IF (judft_was_argument("-info")) THEN
            CLOSE (oUnit)
            OPEN (oUnit, status='SCRATCH')
         ELSE
            inquire (file="out.history", exist=l_exist)
            inquire (file="out", exist=l_found)
            if (l_exist .and. l_found) THEN
               open (666, file="out.history", access="append", status="old")
               open (667, file="out", status="old")
               do
                  read (667, '(a)', end=999) line
                  write (666, '(a)') line
               end do
999            close (667)
               close (666)
            end if
            IF (.NOT. judft_was_argument("-no_out")) &
               OPEN (oUnit, file='out', form='formatted', status='unknown')
         END IF
         CALL writeOutHeader()
         !this should be removed, it deletes output of old inf file
         OPEN (16, status='SCRATCH')
      END IF

      ALLOCATE (t_xcpot_inbuild::xcpot)
      !Only PE==0 reads the fi%input and does basic postprocessing
      IF (fmpi%irank .EQ. 0) THEN
         CALL fleurinput_read_xml(outxmlFileID, filename_add_loc, cell=fi%cell, sym=fi%sym, atoms=fi%atoms, input=fi%input, noco=fi%noco, vacuum=fi%vacuum, field=fi%field, &
                                  sliceplot=fi%sliceplot, banddos=fi%banddos, mpinp=fi%mpinp, hybinp=fi%hybinp, coreSpecInput=fi%coreSpecInput, &
                                  wann=wann, xcpot=xcpot, forcetheo_data=forcetheo_data, kpts=fi%kpts, kptsSelection=kptsSelection, kptsArray=kptsArray, &
                                  enparaXML=enparaXML, gfinp=fi%gfinp, hub1inp=fi%hub1inp, juPhon=fi%juPhon)
         CALL fleurinput_postprocess(fi%cell, fi%sym, fi%atoms, fi%input, fi%noco, fi%vacuum, &
                                     fi%banddos,   Xcpot, fi%kpts, fi%gfinp)
      END IF
      !Distribute fi%input to all PE
      CALL fleurinput_mpi_bc(fi%cell, fi%sym, fi%atoms, fi%input, fi%noco, fi%vacuum, fi%field, &
                             fi%sliceplot, fi%banddos, fi%mpinp, fi%hybinp,   fi%coreSpecInput, Wann, &
                             Xcpot, Forcetheo_data, fi%kpts, Enparaxml, fi%gfinp, fi%hub1inp, fmpi%Mpi_comm, fi%juPhon)
      !Remaining init is done using all PE
      call make_xcpot(fmpi, xcpot, fi%atoms, fi%input)
      IF (fmpi%irank .EQ. 0) THEN
         IF(xcpot%is_hybrid().AND.fi%sym%invs) THEN
            CALL juDFT_warn("Hybrid functionals with inversion symmetric unit cells are disabled at the moment.")
         END IF
      END IF
      CALL nococonv%init(fi%noco)
      CALL nococonv%init_ss(fi%noco, fi%atoms)
      !CALL ylmnorm_init(MAX(fi%atoms%lmaxd, 2*fi%hybinp%lexp))
      CALL gaunt_init(fi%atoms%lmaxd + 1)
      CALL enpara%init_enpara(fi%atoms, fi%input%jspins, fi%input%film, enparaXML)
      CALL make_sphhar(fmpi%irank == 0, fi%atoms, sphhar, fi%sym, fi%cell)
      ! Store structure data (has to be performed before calling make_stars)
      CALL storeStructureIfNew(fi%input, stars, fi%atoms, fi%cell, fi%vacuum,  fi%sym, fmpi, sphhar, fi%noco)
      CALL make_stars(stars, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi)
      CALL make_forcetheo(forcetheo_data, fi%cell, fi%sym, fi%atoms, forcetheo)
      CALL lapw_dim(fi%kpts, fi%cell, fi%input, fi%noco, nococonv,   forcetheo, fi%atoms, nbasfcn)
      CALL fi%input%init(fi%noco, fi%hybinp%l_hybrid,fi%sym%invs,fi%atoms%n_denmat,fi%atoms%n_hia,lapw_dim_nbasfcn)
      CALL fi%hybinp%init(fi%atoms, fi%cell, fi%input,   fi%sym, xcpot)
      l_timeReversalCheck = .FALSE.
      IF(.NOT.fi%banddos%band.AND..NOT.fi%banddos%dos) THEN
         IF(fi%noco%l_soc.OR.fi%noco%l_ss) l_timeReversalCheck = .TRUE.
      END IF
      CALL fi%kpts%init(fi%sym, fi%input%film, fi%hybinp%l_hybrid .or. fi%input%l_rdmft, l_timeReversalCheck)
      CALL fi%kpts%initTetra(fi%input, fi%cell, fi%sym, fi%noco%l_soc .OR. fi%noco%l_ss)
      IF (fmpi%irank == 0) CALL fi%gfinp%init(fi%atoms, fi%sym, fi%noco, fi%cell, fi%input)
      CALL fi%gfinp%mpi_bc(fmpi%mpi_comm) !THis has to be rebroadcasted because there could be new gf elements after init_gfinp
      CALL convn(fmpi%irank == 0, fi%atoms, stars)
      IF (fmpi%irank == 0) CALL e_field(fi%atoms, stars, fi%sym, fi%vacuum, fi%cell, fi%input, fi%field%efield)
      IF (fmpi%isize > 1) CALL fi%field%mpi_bc(fmpi%mpi_comm, 0)

      !At some point this should be enabled for fi%noco as well
      IF (.NOT. fi%noco%l_noco) &
         CALL transform_by_moving_atoms(fmpi, stars,fi%atoms, fi%vacuum, fi%cell, fi%sym, sphhar, fi%input,   fi%noco, nococonv)

      !
      !--> determine more dimensions
      !

      IF (fmpi%irank .EQ. 0) THEN
         CALL writeOutParameters(fmpi, fi%input, fi%sym, stars, fi%atoms, fi%vacuum, fi%kpts, &
                                   fi%hybinp, fi%cell, fi%banddos, fi%sliceplot, xcpot, &
                                 fi%noco, enpara, sphhar)
         CALL fleur_info(fi%kpts)
         CALL deleteDensities()
      END IF

      !Finalize the fmpi setup
      CALL setupMPI(fi%kpts%nkpt, fi%input%neig, nbasfcn, fmpi)

      !Collect some usage info
      CALL add_usage_data("A-Types", fi%atoms%ntype)
      CALL add_usage_data("fi%atoms", fi%atoms%nat)
      CALL add_usage_data("Real", fi%input%l_real)
      CALL add_usage_data("Spins", fi%input%jspins)
      CALL add_usage_data("Noco", fi%noco%l_noco)
      CALL add_usage_data("SOC", fi%noco%l_soc)
      CALL add_usage_data("SpinSpiral", fi%noco%l_ss)
      CALL add_usage_data("PlaneWaves", lapw_dim_nvd)
      CALL add_usage_data("LOs", fi%atoms%nlotot)
      CALL add_usage_data("nkpt", fi%kpts%nkpt)

#ifdef CPP_GPU
      CALL add_usage_data("gpu_per_node", 1)
#else
      CALL add_usage_data("gpu_per_node", 0)
#endif

      CALL results%init(fi%input, fi%atoms, fi%kpts, fi%noco)

      IF (fmpi%irank .EQ. 0) THEN
         IF (fi%input%gw .NE. 0) CALL mixing_history_reset(fmpi)
         CALL setStartingDensity(fi%noco%l_noco)
      END IF

      if(fi%hybinp%l_hybrid) call load_hybrid_data(fi, fmpi, hybdat, mpdata)

      !new check mode will only run the init-part of FLEUR
      IF (judft_was_argument("-check")) THEN  
         call test_performance()
         CALL judft_end("Check-mode done", fmpi%irank)
      endif   
#ifdef CPP_MPI
      CALL MPI_BARRIER(fmpi%mpi_comm, ierr(1))
#endif
   CONTAINS
      SUBROUTINE init_wannier()
         ! Initializations for Wannier functions (start)
         IF (fmpi%irank .EQ. 0) THEN
            wann%l_gwf = wann%l_ms .OR. wann%l_sgwf .OR. wann%l_socgwf

            IF (wann%l_gwf) THEN
               WRITE (*, *) 'running HDWF-extension of FLEUR code'
               WRITE (*, *) 'with l_sgwf =', wann%l_sgwf, ' and l_socgwf =', wann%l_socgwf

               IF (wann%l_socgwf .AND. .NOT. fi%noco%l_soc) THEN
                  CALL juDFT_error("set l_soc=T if l_socgwf=T", calledby="fleur_init")
               END IF

               IF ((wann%l_ms .OR. wann%l_sgwf) .AND. .NOT. (fi%noco%l_noco .AND. fi%noco%l_ss)) THEN
                  CALL juDFT_error("set l_noco=l_ss=T for l_sgwf.or.l_ms", calledby="fleur_init")
               END IF

               IF ((wann%l_ms .OR. wann%l_sgwf) .AND. wann%l_socgwf) THEN
                  CALL juDFT_error("(l_ms.or.l_sgwf).and.l_socgwf", calledby="fleur_init")
               END IF

               INQUIRE (FILE=wann%param_file, EXIST=l_exist)
               IF (.NOT. l_exist) THEN
                  CALL juDFT_error("where is param_file"//TRIM(wann%param_file)//"?", calledby="fleur_init")
               END IF
               OPEN (113, file=wann%param_file, status='old')
               READ (113, *) wann%nparampts, wann%scale_param
               CLOSE (113)
            ELSE
               wann%nparampts = 1
               wann%scale_param = 1.0
            END IF
         END IF

         ALLOCATE (wann%param_vec(3, wann%nparampts))
         ALLOCATE (wann%param_alpha(fi%atoms%ntype, wann%nparampts))

         IF (fmpi%irank .EQ. 0) THEN
            IF (wann%l_gwf) THEN
               OPEN (113, file=wann%param_file, status='old')
               READ (113, *)!header
               WRITE (oUnit, *) 'parameter points for HDWFs generation:'
               IF (wann%l_sgwf .OR. wann%l_ms) THEN
                  WRITE (oUnit, *) '      q1       ', '      q2       ', '      q3'
               ELSE IF (wann%l_socgwf) THEN
                  WRITE (oUnit, *) '      --       ', '     phi       ', '    theta'
               END IF

               DO pc = 1, wann%nparampts
                  READ (113, '(3(f14.10,1x))') wann%param_vec(1, pc), wann%param_vec(2, pc), wann%param_vec(3, pc)
                  wann%param_vec(:, pc) = wann%param_vec(:, pc)/wann%scale_param
                  WRITE (oUnit, '(3(f14.10,1x))') wann%param_vec(1, pc), wann%param_vec(2, pc), wann%param_vec(3, pc)
                  IF (wann%l_sgwf .OR. wann%l_ms) THEN
                     iAtom = 1
                     DO iType = 1, fi%atoms%ntype
                        phi_add = tpi_const*(wann%param_vec(1, pc)*fi%atoms%taual(1, iAtom) + &
                                             wann%param_vec(2, pc)*fi%atoms%taual(2, iAtom) + &
                                             wann%param_vec(3, pc)*fi%atoms%taual(3, iAtom))
                        wann%param_alpha(iType, pc) = nococonv%alph(iType) + phi_add
                        iAtom = iAtom + fi%atoms%neq(iType)
                     END DO
                  END IF
               END DO

               IF (ANY(wann%param_vec(1, :) .NE. wann%param_vec(1, 1))) wann%l_dim(1) = .TRUE.
               IF (ANY(wann%param_vec(2, :) .NE. wann%param_vec(2, 1))) wann%l_dim(2) = .TRUE.
               IF (ANY(wann%param_vec(3, :) .NE. wann%param_vec(3, 1))) wann%l_dim(3) = .TRUE.

               CLOSE (113)

               IF (wann%l_dim(1) .AND. wann%l_socgwf) THEN
                  CALL juDFT_error("do not specify 1st component if l_socgwf", calledby="fleur_init")
               END IF
            END IF!(wann%l_gwf)
         END IF!(fmpi%irank.EQ.0)

#ifdef CPP_MPI
         CALL MPI_BCAST(wann%param_vec, 3*wann%nparampts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr(1))
         CALL MPI_BCAST(wann%param_alpha, fi%atoms%ntype*wann%nparampts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr(1))
         CALL MPI_BCAST(wann%l_dim, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr(1))
#endif

         ! Initializations for Wannier functions (end)
      END SUBROUTINE init_wannier
   END SUBROUTINE fleur_init
END MODULE m_fleur_init
