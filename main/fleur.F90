!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur
   IMPLICIT NONE
CONTAINS
   SUBROUTINE fleur_execute(fmpi, fi, sphhar, stars, nococonv, forcetheo, enpara, results, &
                            xcpot, wann, hybdat, mpdata)

      !     ***************************************************************
      !
      !     based on flapw7 (c.l.fu, m.weinert, e.wimmer):
      !     full potential linearized augmented plane wave method for thin
      !     films and superlattices (version 7 ---- general symmetry)
      !     symmetry part       ---  e.wimmer
      !     potential generator ---  c.l.fu,r.podloucky
      !     matrix elements     ---  m.weinert
      !     charge density      ---  c.l.fu
      !                                c.l.fu        1987
      !     2nd variation diagon.  --- r.-q. wu      1992
      !     forces a la Yu et al   --- r.podloucky   1995
      !     full relativistic core --- a.shick       1996
      !     broyden mixing         --- r.pentcheva   1996
      !     gga (pw91, pbe)        --- t.asada       1997
      !     local orbitals         --- p.kurz        1997
      !     automatic symmetry     --- w.hofer       1997
      !     core tails & start     --- r.abt         1998
      !     spin orbit coupling    --- a.shick,x.nie 1998
      !     non-colinear magnet.   --- p.kurz        1999
      !     one-dimensional        --- y.mokrousov   2002
      !     exchange parameters    --- m.lezaic      2004
      !
      !                       g.bihlmayer, s.bluegel 1999
      !     ***************************************************************
      !----------------------------------------
      ! this routine is the main PROGRAM

      USE m_types
      USE m_types_forcetheo_extended
      USE m_constants
      USE m_optional
      USE m_cdn_io
      USE m_mixing_history
      USE m_qfix
      USE m_vgen
      USE m_vgen_coulomb
      USE m_writexcstuff
      USE m_eigen
      USE m_eigenso
      USE m_fermie
      USE m_cdngen
      USE m_totale
      USE m_potdis
      USE m_mix
      USE m_xmlOutput
      USE m_juDFT_time
      USE m_calc_hybrid
      USE m_rdmft
      USE m_io_hybrid
      USE m_wann_optional
      USE m_wannier
      USE m_bs_comfort
      USE m_dwigner
      USE m_ylm
      USE m_metagga
      USE m_plot
      USE m_usetup
      USE m_hubbard1_setup
      USE m_writeCFOutput
      USE m_mpi_bc_tool
      USE m_eig66_io
      USE m_chase_diag
      USE m_writeBasis
      USE m_RelaxSpinAxisMagn
      USE m_dfpt

!$    USE omp_lib
      IMPLICIT NONE

      TYPE(t_mpi), INTENT(INOUT)      :: fmpi
      type(t_fleurinput), intent(in) :: fi
      class(t_xcpot), intent(in)     :: xcpot
      TYPE(t_sphhar), INTENT(IN)      :: sphhar
      TYPE(t_stars), INTENT(IN)       :: stars
      TYPE(t_nococonv), intent(inout) :: nococonv
      type(t_results), intent(inout) :: results
      type(t_wann), intent(inout)    :: wann

      CLASS(t_forcetheo), INTENT(INOUT)::forcetheo
      TYPE(t_enpara), INTENT(INOUT)    :: enpara
      TYPE(t_hybdat), intent(inout)    :: hybdat

      TYPE(t_input) :: input_soc !same as fi%input with neig=2*neig !should be refactored out

      TYPE(t_field)                    :: field2
      TYPE(t_mpdata), intent(inout)    :: mpdata

      TYPE(t_potden)                  :: vTot, vx, vCoul, vxc, exc
      TYPE(t_potden)                  :: inDen, outDen, EnergyDen, sliceDen

      TYPE(t_hub1data)                :: hub1data
      TYPE(t_greensf), ALLOCATABLE    :: greensFunction(:)

      ! local scalars
      INTEGER :: eig_id, archiveType, num_threads
      INTEGER :: iter, iterHF, i, n, i_gf
      INTEGER :: wannierspin
      LOGICAL :: l_opti, l_cont, l_qfix, l_real, l_olap, l_error, l_dummy
      LOGICAL :: l_forceTheorem, l_lastIter
      REAL    :: fix, sfscale, rdummy, tempDistance
      REAL    :: mmpmatDistancePrev,occDistancePrev

#ifdef CPP_MPI
      INTEGER :: ierr
#endif
      REAL, ALLOCATABLE :: flh(:, :), flh2(:, :)
      COMPLEX, ALLOCATABLE :: flm(:, :)
      INTEGER, ALLOCATABLE :: nvfull(:, :), GbasVec_eig(:, :, :, :)

      IF ((fi%input%preconditioning_param /= 0) .AND. fi%oneD%odi%d1) THEN
         CALL juDFT_error('Currently no preconditioner for 1D calculations', calledby='fleur')
      END IF

      CALL optional(fmpi, fi%atoms, sphhar, fi%vacuum, &
                    stars, fi%input, fi%sym, fi%cell, fi%sliceplot, xcpot, fi%noco, fi%oneD)

      IF (fi%input%l_wann .AND. (fmpi%irank == 0) .AND. (.NOT. wann%l_bs_comf)) THEN
!       IF(fmpi%isize.NE.1) CALL juDFT_error('No Wannier+MPI at the moment',calledby = 'fleur')
         CALL wann_optional(fmpi, fi%input, fi%kpts, fi%atoms, fi%sym, fi%cell, fi%oneD, fi%noco, wann)
      END IF

      iter = 0
      iterHF = 0
      l_cont = (iter < fi%input%itmax)

      !Read in last Hubbard 1 distances
      l_error = .TRUE.
      IF(fi%atoms%n_hia>0 .AND. fmpi%irank.EQ.0) CALL readPrevmmpDistances(mmpmatDistancePrev,occDistancePrev,l_error)
      CALL hub1data%init(fi%atoms, fi%input, fi%hub1inp, fmpi, mmpmatDistancePrev, occDistancePrev, l_error)
      CALL hub1data%mpi_bc(fmpi%mpi_comm)
      IF(fi%atoms%n_hia>0 .AND. .NOT.l_error .AND. .NOT.fi%hub1inp%l_forceHIAiteration) THEN
         !Set the current HIA distance to the read in value
         !Prevents too many HIA iterations after restart
         results%last_mmpmatDistance = mmpmatDistancePrev
         results%last_occDistance = occDistancePrev
      ENDIF
      CALL mpi_bc(results%last_mmpmatDistance,0,fmpi%mpi_comm)
      CALL mpi_bc(results%last_occDistance,0,fmpi%mpi_comm)

      IF (fmpi%irank .EQ. 0) CALL openXMLElementNoAttributes('scfLoop')

      ! Initialize and load inDen density (start)

      !Warning on strange choice of switches before starting density is generated.
      IF (fi%input%l_onlyMtStDen .AND. .NOT. any(fi%noco%l_unrestrictMT)) THEN
         CALL juDFT_warn("l_onlyMtStDen='T' and l_mtNocoPot='F' makes no sense.", calledby='types_input')
      END IF

      CALL inDen%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN)

      archiveType = CDN_ARCHIVE_TYPE_CDN1_const
      IF (fi%noco%l_noco) archiveType = CDN_ARCHIVE_TYPE_NOCO_const
      IF (any(fi%noco%l_unrestrictMT)) archiveType = CDN_ARCHIVE_TYPE_FFN_const
      IF (fmpi%irank .EQ. 0) CALL readDensity(stars, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, fi%input, fi%sym, fi%oneD, archiveType, CDN_INPUT_DEN_const, &
                                              0, results%ef, results%last_distance, l_qfix, inDen)
      !IF (fi%noco%l_alignMT .AND. fmpi%irank .EQ. 0) THEN
      !   CALL initRelax(fi%noco, nococonv, fi%atoms, fi%input, fi%vacuum, sphhar, stars, fi%sym, fi%oneD, fi%cell, inDen)
         !CALL doRelax(fi%vacuum, sphhar, stars, fi%sym, fi%oneD, fi%cell, fi%noco, nococonv, fi%input, fi%atoms, inDen)
      !END IF
      CALL timestart("Qfix")
      CALL qfix(fmpi, stars, fi%atoms, fi%sym, fi%vacuum, sphhar, fi%input, fi%cell, fi%oneD, inDen, fi%noco%l_noco, .FALSE., .FALSE., .FALSE., fix)
      CALL timestop("Qfix")
      IF (fmpi%irank .EQ. 0) THEN
         CALL writeDensity(stars, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, fi%input, fi%sym, fi%oneD, archiveType, CDN_INPUT_DEN_const, &
                           0, -1.0, results%ef, results%last_mmpmatDistance, results%last_occDistance, .FALSE., inDen)
      END IF

      IF (any(fi%noco%l_alignMT)) CALL toLocalSpinFrame(fmpi,fi%vacuum, sphhar, stars, fi%sym, fi%oneD, fi%cell, fi%noco, nococonv, fi%input, fi%atoms, .true.,inDen,.true.)
      ! Initialize and load inDen density (end)

      ! Initialize potentials (start)
      CALL vTot%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL vCoul%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTCOUL)
      CALL vx%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL vxc%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL exc%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT)
      ! Initialize potentials (end)

      ! Initialize Green's function (start)
      ALLOCATE (greensFunction(MAX(1, fi%gfinp%n)))
      IF (fi%gfinp%n > 0) THEN
         DO i_gf = 1, fi%gfinp%n
            CALL greensFunction(i_gf)%init(fi%gfinp%elem(i_gf), fi%gfinp, fi%atoms, fi%input, nkpt=fi%kpts%nkpt)
         ENDDO
      ENDIF
      ! Initialize Green's function (end)

      ! Open/allocate eigenvector storage (start)
      if (fi%noco%l_soc .and. fi%input%l_wann) then
       !! Weed up and down spinor components for SOC MLWFs.
       !! When jspins=1 Fleur usually writes only the up-spinor into the eig-file.
       !! Make sure we always get up and down spinors when SOC=true.
         wannierspin = 2
      else
         wannierspin = fi%input%jspins
      endif
      l_olap = fi%hybinp%l_hybrid .OR. fi%input%l_rdmft
      eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, wannierspin, &
                        fi%noco%l_noco,.NOT. fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%INPUT%eig66(1), l_olap, fmpi%n_size)
      hybdat%eig_id = eig_id
!Rotate cdn to local frame if specified.

#ifdef CPP_CHASE
      CALL init_chase(fmpi, fi%input, fi%atoms, fi%kpts, fi%noco, l_real)
#endif
      ! Open/allocate eigenvector storage (end)

      ! Perform some pre-scf-loop checks (start)
      IF (fi%input%gw .EQ. 2) THEN
         IF (ABS(results%last_distance).GT.fi%input%mindistance) THEN
            CALL juDFT_warn('Performing spex="2" step without a selfconsistent density!')
         END IF
      END IF
      ! Perform some pre-scf-loop checks (end)

      l_lastIter = .FALSE.
      scfloop: DO WHILE (l_cont)
         iter = iter + 1
         l_lastIter = l_lastIter.OR.(iter.EQ.fi%input%itmax)
         hub1data%overallIteration = hub1data%overallIteration + 1

         IF (fmpi%irank .EQ. 0) CALL openXMLElementFormPoly('iteration', (/'numberForCurrentRun', 'overallNumber      '/), &
                                                            (/iter, inden%iter/), RESHAPE((/19, 13, 5, 5/), (/2, 2/)))

!!$       !+t3e
!!$       IF (fi%input%alpha.LT.10.0) THEN
!!$
!!$          IF (iter.GT.1) THEN
!!$             fi%input%alpha = fi%input%alpha - NINT(fi%input%alpha)
!!$          END IF

         !CALL resetIterationDependentTimers()
         CALL timestart("Iteration")
         IF (fmpi%irank .EQ. 0) THEN
            WRITE (oUnit, FMT=8100) iter
8100        FORMAT(/, 10x, '   iter=  ', i5)
         ENDIF !fmpi%irank.eq.0

#ifdef CPP_CHASE
         CALL chase_distance(results%last_distance)
#endif

         CALL inDen%distribute(fmpi%mpi_comm)
         call nococonv%mpi_bc(fmpi%mpi_comm)

!Plot inden if wanted
         IF (fi%sliceplot%iplot .NE. 0) THEN
            IF (.NOT.fi%sliceplot%slice) THEN
               CALL makeplots(stars, fi%atoms, sphhar, fi%vacuum, fi%input, fmpi, fi%oneD, fi%sym, fi%cell, &
                              fi%noco, nococonv, inDen, PLOT_INPDEN, fi%sliceplot)
            ELSE
               CALL sliceDen%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN)
               IF (fmpi%irank .EQ. 0) CALL readDensity(stars, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, &
                                                       fi%input, fi%sym, fi%oneD,CDN_ARCHIVE_TYPE_CDN_const, &
                                                       CDN_INPUT_DEN_const, 0, rdummy, tempDistance, l_dummy, sliceDen, 'cdn_slice')
               CALL sliceden%distribute(fmpi%mpi_comm)
               call nococonv%mpi_bc(fmpi%mpi_comm)
               CALL makeplots(stars, fi%atoms, sphhar, fi%vacuum, fi%input, fmpi, fi%oneD, fi%sym, fi%cell, &
                              fi%noco, nococonv, sliceDen, PLOT_INPDEN, fi%sliceplot)
            END IF
            IF ((fmpi%irank .EQ. 0) .AND. (fi%sliceplot%iplot .EQ. 2)) THEN
               CALL juDFT_end("Stopped self consistency loop after plots have been generated.")
            END IF
         END IF

         !HF
         IF (fi%hybinp%l_hybrid) THEN
            SELECT TYPE (xcpot)
            TYPE IS (t_xcpot_inbuild)
               CALL calc_hybrid(fi, mpdata, hybdat, fmpi, nococonv, stars, enpara, &
                                results, xcpot, vTot, iter, iterHF)
            END SELECT
#ifdef CPP_MPI
            call MPI_Barrier(fmpi%mpi_comm, ierr)
#endif
            IF (hybdat%l_calhf) THEN
               call mixing_history_reset(fmpi)
               iter = 0
            END IF
         ENDIF

!!$             DO pc = 1, wann%nparampts
!!$                !---> gwf
!!$                IF (wann%l_sgwf.OR.wann%l_ms) THEN
!!$                   fi%noco%qss(:) = wann%param_vec(:,pc)
!!$                   fi%noco%alph(:) = wann%param_alpha(:,pc)
!!$                ELSE IF (wann%l_socgwf) THEN
!!$                   IF(wann%l_dim(2)) fi%noco%phi   = tpi_const * wann%param_vec(2,pc)
!!$                   IF(wann%l_dim(3)) fi%noco%theta = tpi_const * wann%param_vec(3,pc)
!!$                END IF
         !---< gwf

         IF (any(fi%noco%l_unrestrictMT) .AND. fi%noco%l_scaleMag) THEN
            sfscale = fi%noco%mag_scale
            CALL inDen%SpinsToChargeAndMagnetisation()
            inDen%mt(:, 0:, :, 2:4) = sfscale*inDen%mt(:, 0:, :, 2:4)
            inDen%pw(:, 2:3) = sfscale*inDen%pw(:, 2:3)
            inDen%vacz(:, :, 2:4) = sfscale*inDen%vacz(:, :, 2:4)
            inDen%vacxy(:, :, :, 2:3) = sfscale*inDen%vacxy(:, :, :, 2:3)
            CALL inDen%ChargeAndMagnetisationToSpins()
         END IF

         CALL timestart("generation of potential")
         CALL vgen(hybdat, fi%field, fi%input, xcpot, fi%atoms, sphhar, stars, fi%vacuum, fi%sym, &
                   fi%cell, fi%oneD, fi%sliceplot, fmpi, results, fi%noco, nococonv, EnergyDen, inDen, vTot, vx, vCoul, vxc, exc)
         CALL timestop("generation of potential")

         IF (any(fi%noco%l_unrestrictMT) .AND. fi%noco%l_scaleMag) THEN
            CALL inDen%SpinsToChargeAndMagnetisation()
            inDen%mt(:, 0:, :, 2:4) = inDen%mt(:, 0:, :, 2:4)/sfscale
            inDen%pw(:, 2:3) = inDen%pw(:, 2:3)/sfscale
            inDen%vacz(:, :, 2:4) = inDen%vacz(:, :, 2:4)/sfscale
            inDen%vacxy(:, :, :, 2:3) = inDen%vacxy(:, :, :, 2:3)/sfscale
            CALL inDen%ChargeAndMagnetisationToSpins()
         END IF


         IF (hub1data%l_runthisiter .AND. fi%atoms%n_hia > 0) THEN
            DO i_gf = 1, fi%gfinp%n
               CALL greensFunction(i_gf)%mpi_bc(fmpi%mpi_comm)
            ENDDO
            IF (ALL(greensFunction(fi%gfinp%hiaElem)%l_calc)) THEN
               hub1data%iter = hub1data%iter + 1
               CALL hubbard1_setup(fi%atoms, fi%cell, fi%gfinp, fi%hub1inp, fi%input, fmpi, fi%noco, fi%kpts, nococonv, vTot, &
                                   greensFunction(fi%gfinp%hiaElem), hub1data, results, inDen)
            ELSE
               IF (fmpi%irank .EQ. 0) WRITE (*, *) 'Not all Greens Functions available: Running additional iteration'
               hub1data%l_runthisiter = .FALSE. !To prevent problems in mixing later on
            ENDIF
         ENDIF

         IF (fi%atoms%n_u + fi%atoms%n_hia > 0) THEN
            CALL u_setup(fi%atoms, fi%input, fi%noco, fmpi, hub1data, inDen, vTot, results)
         END IF

#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif
         CALL forcetheo%start(vtot, fmpi%irank == 0)
         forcetheoloop: DO WHILE (forcetheo%next_job(fmpi,l_lastIter, fi%atoms, fi%noco, nococonv))

            CALL timestart("gen. of hamil. and diag. (total)")
            CALL timestart("eigen")
            CALL timestart("Updating energy parameters")
            CALL enpara%update(fmpi%mpi_comm, fi%atoms, fi%vacuum, fi%input, vToT, hub1data)
            CALL timestop("Updating energy parameters")
            IF (.not. fi%input%eig66(1)) THEN
                IF (fi%juPhon%l_dfpt) THEN
                    CALL eigen(fi, fmpi, stars, sphhar, xcpot, &
                               enpara, nococonv, mpdata, hybdat, &
                               iter, eig_id, results, inDen, vToT, vx, hub1data, nvfull, GbasVec_eig)
                ELSE
                    CALL eigen(fi, fmpi, stars, sphhar, xcpot, &
                               enpara, nococonv, mpdata, hybdat, &
                               iter, eig_id, results, inDen, vToT, vx, hub1data)
                END IF
            ENDIF
!!$          eig_idList(pc) = eig_id
            CALL timestop("eigen")

            ! add all contributions to total energy
#ifdef CPP_MPI
            ! send all result of local total energies to the r
            IF (fi%hybinp%l_hybrid .AND. hybdat%l_calhf) THEN
               call fmpi%set_root_comm()
               if (fmpi%n_rank == 0) then
                  IF (fmpi%irank == 0) THEN
                     CALL MPI_Reduce(MPI_IN_PLACE, results%te_hfex%valence, 1, MPI_REAL8, MPI_SUM, 0, fmpi%root_comm, ierr)
                  ELSE
                     CALL MPI_Reduce(results%te_hfex%valence, MPI_IN_PLACE, 1, MPI_REAL8, MPI_SUM, 0, fmpi%root_comm, ierr)
                  END IF
               endif
            END IF
#endif

            ! WRITE(oUnit,fmt='(A)') 'Starting 2nd variation ...'
            IF (fi%noco%l_soc .AND. .NOT. fi%noco%l_noco .AND. .NOT. fi%INPUT%eig66(1)) &
               CALL eigenso(eig_id, fmpi, stars, sphhar, nococonv, vTot, enpara, results, fi%hub1inp, hub1data,fi)
            CALL timestop("gen. of hamil. and diag. (total)")

#ifdef CPP_MPI
            CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

            ! fermi level and occupancies
            input_soc = fi%input
            IF (fi%noco%l_soc .AND. (.NOT. fi%noco%l_noco)) then
               input_soc = fi%input
               input_soc%neig = 2*fi%input%neig
            ENDIF

            IF (fi%input%gw .GT. 0) THEN
               IF (fmpi%irank .EQ. 0) THEN
                  CALL writeBasis(input_soc, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, fi%cell, enpara, hub1data, vTot, vCoul, vx, fmpi, &
                                  results, eig_id, fi%oneD, sphhar, stars, fi%vacuum)
               END IF
               IF (fi%input%gw .EQ. 2) THEN
                  CALL juDFT_end("SPEX data written. Fleur ends.", fmpi%irank)
               END IF
            END IF

            !IF ((fmpi%irank.EQ.0)) THEN
            CALL timestart("determination of fermi energy")

            IF (fi%noco%l_soc .AND. (.NOT. fi%noco%l_noco)) THEN
               input_soc%zelec = fi%input%zelec*2
               CALL fermie(eig_id, fmpi, fi%kpts, input_soc, fi%noco, enpara%epara_min, fi%cell, results)
               results%seigscv = results%seigscv/2
               results%ts = results%ts/2
            ELSE
               CALL fermie(eig_id, fmpi, fi%kpts, fi%input, fi%noco, enpara%epara_min, fi%cell, results)
            ENDIF
            CALL timestop("determination of fermi energy")

!!$          !+Wannier
!!$          IF(wann%l_bs_comf)THEN
!!$             IF(pc.EQ.1) THEN
!!$                OPEN(777,file='out_eig.1')
!!$                OPEN(778,file='out_eig.2')
!!$                OPEN(779,file='out_eig.1_diag')
!!$                OPEN(780,file='out_eig.2_diag')
!!$             END IF
!!$
!!$             CALL bs_comfort(eig_id,fi%input,fi%noco,fi%kpts%nkpt,pc)
!!$
!!$             IF(pc.EQ.wann%nparampts)THEN
!!$                CLOSE(777)
!!$                CLOSE(778)
!!$                CLOSE(779)
!!$                CLOSE(780)
!!$             END IF
!!$          END IF
!!$          !-Wannier

            !ENDIF
#ifdef CPP_MPI
            CALL MPI_BCAST(results%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            CALL MPI_BCAST(results%w_iks, SIZE(results%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif

            !!!juPhon:
            ! Implement here a call to a new jpMain. In it, generate all necessary
            ! quantities only from the density (potentials), transform them into an
            ! unsymmetrized form and basically cover all steps that were formerly in
            ! jpInit. Then do the whole juPhon shtick until we have phonons. This
            ! eliminates the need for complicated juggling with the eig file and
            ! unoccupied states.
            !
            ! Also add safety. We do not want to even enter juPhon if ctail is on,
            ! jspins is 2, natoms is 2 or more, noco is on, libxc is off etc. etc.
            !
            ! Forbidden switches:
            ! symor = .FALSE.
            !
            ! Those restrictions will gradually be lifted.

            IF (fi%juPhon%l_dfpt) THEN
                CALL dfpt(fi%juPhon, fi%sym, fi%input, fi%atoms, sphhar, stars, fi%cell, fi%noco, nococonv, &
                        & fi%kpts, fmpi, results, enpara, inDen, vTot, eig_id, nvfull, GbasVec_eig)
            END IF

            !!!juPhon

            IF (forcetheo%eval(eig_id, fi%atoms, fi%kpts, fi%sym, fi%cell, fi%noco, nococonv, input_soc, fmpi, fi%oneD, enpara, vToT, results)) THEN
               CYCLE forcetheoloop
            ENDIF

            !+Wannier functions
            IF ((fi%input%l_wann) .AND. (.NOT. wann%l_bs_comf)) THEN
               CALL wannier(fmpi, input_soc, fi%kpts, fi%sym, fi%atoms, stars, fi%vacuum, sphhar, fi%oneD, &
                            wann, fi%noco, nococonv, fi%cell, enpara, fi%banddos, fi%sliceplot, vTot, results, &
                            (/eig_id/), (fi%sym%invs) .AND. (.NOT. fi%noco%l_soc) .AND. (.NOT. fi%noco%l_noco), fi%kpts%nkpt)
            END IF
            !-Wannier

            !Check if the greensFunction have to be calculated
            IF (fi%gfinp%n > 0) THEN
               DO i_gf = 1, fi%gfinp%n
                  !Either the set distance has been reached (or is negative)
                  greensFunction(i_gf)%l_calc = (results%last_distance >= 0.0 .AND. &
                                                 results%last_distance < fi%gfinp%minCalcDistance) &
                                                .OR. fi%gfinp%minCalcDistance < 0.0 & !No minCalcDistance distance set
                                                .OR. iter == fi%input%itmax !Maximum iteration  reached
                  !or we are in the first iteration for Hubbard 1
                  IF (fi%atoms%n_hia > 0) THEN
                     greensFunction(i_gf)%l_calc = greensFunction(i_gf)%l_calc .OR. (iter == 1 .AND. (hub1data%iter == 0 &
                                                                                                      .AND. ALL(ABS(vTot%mmpMat(:, :, fi%atoms%n_u + 1:fi%atoms%n_u + fi%atoms%n_hia, :)) .LT. 1e-12)))
                  ENDIF
               ENDDO
            ENDIF

            ! charge density generation
            CALL timestart("generation of new charge density (total)")
            CALL outDen%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN)
            outDen%iter = inDen%iter
            CALL cdngen(eig_id, fmpi, input_soc, fi%banddos, fi%sliceplot, fi%vacuum, &
                        fi%kpts, fi%atoms, sphhar, stars, fi%sym, fi%gfinp, fi%hub1inp, &
                        enpara, fi%cell, fi%noco, nococonv, vTot, results, fi%oneD, fi%corespecinput, &
                        archiveType, xcpot, outDen, EnergyDen, greensFunction, hub1data,vxc,exc)
            !The density matrix for DFT+Hubbard1 only changes in hubbard1_setup and is kept constant otherwise
            outDen%mmpMat(:, :, fi%atoms%n_u + 1:fi%atoms%n_u + fi%atoms%n_hia, :) = inDen%mmpMat(:, :, fi%atoms%n_u + 1:fi%atoms%n_u + fi%atoms%n_hia, :)

            IF (fi%sliceplot%iplot .NE. 0) THEN
               !               CDN including core charge

               CALL makeplots(stars, fi%atoms, sphhar, fi%vacuum, fi%input, fmpi, fi%oneD, fi%sym, &
                              fi%cell, fi%noco, nococonv, outDen, PLOT_OUTDEN_Y_CORE, fi%sliceplot)

               IF ((fi%sliceplot%iplot .NE. 0) .AND. (fmpi%irank .EQ. 0) .AND. (fi%sliceplot%iplot .LT. 64) .AND. (MODULO(fi%sliceplot%iplot, 2) .NE. 1)) THEN
                  CALL juDFT_end("Stopped self consistency loop after plots have been generated.")
               END IF

                END IF

            IF (fi%input%l_rdmft) THEN
               SELECT TYPE (xcpot)
               TYPE IS (t_xcpot_inbuild)
                  CALL rdmft(eig_id, fmpi, fi, enpara, stars, &
                             sphhar, vTot, vCoul, nococonv, xcpot, mpdata, hybdat, &
                             results, archiveType, outDen)
               END SELECT
            END IF

#ifdef CPP_MPI
            CALL MPI_BCAST(enpara%evac, SIZE(enpara%evac), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            CALL MPI_BCAST(enpara%evac0, SIZE(enpara%evac0), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            CALL MPI_BCAST(enpara%el0, SIZE(enpara%el0), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            CALL MPI_BCAST(enpara%ello0, SIZE(enpara%ello0), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)

            !IF (fi%noco%l_noco) THEN
            !   DO n = 1, fi%atoms%ntype
            !      IF (fi%noco%l_relax(n)) THEN
            !         CALL MPI_BCAST(nococonv%alph(n), 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            !         CALL MPI_BCAST(nococonv%beta(n), 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            !      ENDIF
            !   ENDDO
            !   IF (any(fi%noco%l_constrained)) THEN
            !      CALL MPI_BCAST(nococonv%b_con, SIZE(nococonv%b_con), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            !   ENDIF
            !ENDIF
#endif
            CALL timestop("generation of new charge density (total)")


            !CRYSTAL FIELD OUTPUT
            IF(ANY(fi%atoms%l_outputCFpot(:)).OR.ANY(fi%atoms%l_outputCFcdn(:))) THEN
               CALL hub1data%mpi_bc(fmpi%mpi_comm)
               CALL writeCFOutput(fi,stars,hybdat,sphhar,xcpot,EnergyDen,outDen,hub1data,nococonv,enpara,fmpi)
               CALL juDFT_end("Crystal Field Output written",fmpi%irank)
            ENDIF


!!$             !----> output potential and potential difference
!!$             IF (disp) THEN
!!$                reap = .FALSE.
!!$                CALL timestart("generation of potential (total)")
!!$                CALL vgen(fi%hybinp,reap,fi%input,xcpot, fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
!!$                     fi%cell,fi%oneD,fi%sliceplot,fmpi, results,fi%noco,outDen,inDenRot,vTot,vx,vCoul)
!!$                CALL timestop("generation of potential (total)")
!!$
!!$                CALL potdis(stars,fi%vacuum,fi%atoms,sphhar, fi%input,fi%cell,fi%sym)
!!$             END IF

            ! total energy

            !Rotating from local MT frame in global frame for mixing
            CALL toGlobalSpinFrame(fi%noco, nococonv, fi%vacuum, sphhar, stars, fi%sym, fi%oneD, fi%cell, fi%input, fi%atoms, inDen,fmpi)
            CALL toGlobalSpinFrame(fi%noco, nococonv, fi%vacuum, sphhar, stars, fi%sym, fi%oneD, fi%cell, fi%input, fi%atoms, outDen,fmpi,.true.)
            CALL timestart('determination of total energy')
            CALL totale(fmpi, fi%atoms, sphhar, stars, fi%vacuum, fi%sym, fi%input, fi%noco, fi%cell, fi%oneD, &
                        xcpot, hybdat, vTot, vCoul, iter, inDen, results)
            CALL timestop('determination of total energy')
         END DO forcetheoloop

         CALL forcetheo%postprocess()

         CALL enpara%mix(fmpi%mpi_comm, fi%atoms, fi%vacuum, fi%input, vTot%mt(:, 0, :, :), vtot%vacz)
         field2 = fi%field
         ! mix fi%input and output densities
         CALL mix_charge(field2, fmpi, (iter == fi%input%itmax .OR. judft_was_argument("-mix_io")), &
         stars, fi%atoms, sphhar, fi%vacuum, fi%input, &
         fi%sym, fi%cell, fi%noco,nococonv, fi%oneD, archiveType, xcpot, iter, inDen, outDen, results, hub1data%l_runthisiter, fi%sliceplot)

         !Rotating in local MT frame
         CALL toLocalSpinFrame(fmpi,fi%vacuum, sphhar, stars &
         , fi%sym, fi%oneD, fi%cell, fi%noco, nococonv, fi%input, fi%atoms, .true.,inDen,.true.)

         IF (fmpi%irank == 0) THEN
            WRITE (oUnit, FMT=8130) iter
8130        FORMAT(/, 5x, '******* it=', i3, '  is completed********', /,/)
            WRITE (*, *) "Iteration:", iter, " Distance:", results%last_distance
         END IF ! fmpi%irank.EQ.0
         CALL timestop("Iteration")

#ifdef CPP_MPI
         CALL MPI_BCAST(results%last_distance, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
         CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif
         CALL priv_geo_end(fmpi)

         l_cont = .TRUE.
         IF (fi%hybinp%l_hybrid) THEN
            IF (hybdat%l_calhf) THEN
               l_cont = l_cont .AND. (iterHF < fi%input%itmax)
               l_cont = l_cont .AND. (fi%input%mindistance <= results%last_distance)
               CALL check_time_for_next_iteration(iterHF, l_cont)
            ELSE
               l_cont = l_cont .AND. (iter < 50) ! Security stop for non-converging nested PBE calculations
            END IF
            IF (hybdat%l_subvxc) THEN
               results%te_hfex%valence = 0
            END IF
         ELSE IF (fi%atoms%n_hia > 0) THEN
            l_cont = l_cont .AND. (iter < fi%input%itmax) !The SCF cycle reached the maximum iteration
            l_cont = l_cont .AND. ((fi%input%mindistance <= results%last_distance) .OR. fi%input%l_f)
            !If we have converged run hia if the density matrix has not converged
            hub1data%l_runthisiter = .NOT. l_cont .AND. (fi%hub1inp%minoccDistance <= results%last_occdistance &
                                                         .OR. results%last_occdistance <= 0.0 .OR. results%last_mmpMatdistance <= 0.0 &
                                                         .OR. fi%hub1inp%minmatDistance <= results%last_mmpMatdistance)
            !Run after first overall iteration to generate a starting density matrix
            hub1data%l_runthisiter = hub1data%l_runthisiter .OR. (iter == 1 .AND. (hub1data%iter == 0 &
                                                            .AND. ALL(ABS(vTot%mmpMat(:, :, fi%atoms%n_u + 1:fi%atoms%n_u + fi%atoms%n_hia, :)) .LT. 1e-12)))
            hub1data%l_runthisiter = hub1data%l_runthisiter .AND. (iter < fi%input%itmax)
            hub1data%l_runthisiter = hub1data%l_runthisiter .AND. (hub1data%iter < fi%hub1inp%itmax)
            !Prevent that the scf loop terminates
            l_cont = l_cont .OR. hub1data%l_runthisiter
            CALL check_time_for_next_iteration(hub1data%overallIteration, l_cont)
         ELSE
            l_cont = l_cont .AND. (iter < fi%input%itmax)
            ! MetaGGAs need a at least 2 iterations
            l_cont = l_cont .AND. ((fi%input%mindistance <= results%last_distance) .OR. fi%input%l_f &
                                   .OR. (xcpot%exc_is_MetaGGA() .and. iter == 1))
            CALL check_time_for_next_iteration(iter, l_cont)
         END IF

         ! Add extra iteration for force theorem if necessary
         l_forceTheorem = .FALSE.
         SELECT TYPE(forcetheo)
            TYPE IS(t_forcetheo_mae)
               l_forceTheorem = .TRUE.
            TYPE IS(t_forcetheo_dmi)
               l_forceTheorem = .TRUE.
            TYPE IS(t_forcetheo_jij)
               l_forceTheorem = .TRUE.
            TYPE IS(t_forcetheo_ssdisp)
               l_forceTheorem = .TRUE.
         END SELECT

         IF(l_forceTheorem.AND..NOT.l_cont) THEN
            IF(.NOT.l_lastIter) THEN
               l_lastIter = .TRUE.
               l_cont = .TRUE.
            END IF
         ELSE IF(l_forceTheorem.AND.l_lastIter) THEN
            l_cont = .FALSE.
         END IF

         !CALL writeTimesXML()

         IF (fmpi%irank .EQ. 0) THEN
            IF (isCurrentXMLElement("iteration")) CALL closeXMLElement('iteration')
         END IF

         ! Plots of mixed density
         IF ((fi%sliceplot%iplot .NE. 0)) THEN
            ! CDN including core charge
            CALL makeplots(stars, fi%atoms, sphhar, fi%vacuum, fi%input, fmpi, fi%oneD, fi%sym, &
                           fi%cell, fi%noco, nococonv, inDen, PLOT_MIXDEN_Y_CORE, fi%sliceplot)
           !! CDN subtracted by core charge
            !CALL makeplots(fi%sym,stars,fi%vacuum,fi%atoms,sphhar,fi%input,fi%cell,fi%oneD,fi%noco,fi%sliceplot,inDen,PLOT_MIXDEN_N_CORE)
            !CALL makeplots(stars, fi%atoms, sphhar, fi%vacuum, fi%input, fi%oneD, fi%sym, &
            !fi%cell, fi%noco, inDen, PLOT_OUTDEN_N_CORE, fi%sliceplot)
         END IF

         ! Break SCF loop if Plots were generated in ongoing run (iplot=/=0). This needs to happen here, as the mixed density
         ! is the last plottable t_potden to appear in the scf loop and with no mixed density written out (so it is quasi
         ! post-process).

         IF ((fi%sliceplot%iplot .NE. 0) .AND. (fmpi%irank .EQ. 0)) THEN
            CALL juDFT_end("Stopped self consistency loop after plots have been generated.")
         END IF

      END DO scfloop ! DO WHILE (l_cont)

      CALL add_usage_data("Iterations", iter)

      IF (fmpi%irank .EQ. 0) CALL closeXMLElement('scfLoop')

      CALL close_eig(eig_id)
      CALL juDFT_end("all done", fmpi%irank)
   CONTAINS
      SUBROUTINE priv_geo_end(fmpi)
         TYPE(t_mpi), INTENT(IN)::fmpi
         LOGICAL :: l_exist
         !Check if a new fi%input was generated
         INQUIRE (file='inp_new', exist=l_exist)
         IF (l_exist) THEN
            CALL juDFT_end(" GEO new inp created ! ", fmpi%irank)
         END IF
         !check for inp.xml
         INQUIRE (file='inp_new.xml', exist=l_exist)
         IF (.NOT. l_exist) RETURN
         IF (fmpi%irank == 0) THEN
            CALL system('mv inp.xml inp_old.xml')
            CALL system('mv inp_new.xml inp.xml')
            INQUIRE (file='qfix', exist=l_exist)
            IF (l_exist) THEN
               CALL juDFT_end(" GEO new inp created ! ", fmpi%irank)
            END IF
            !check for inp.xml
            INQUIRE (file='inp_new.xml', exist=l_exist)
            IF (.NOT. l_exist) RETURN
            IF (fmpi%irank == 0) THEN
               CALL system('mv inp.xml inp_old.xml')
               CALL system('mv inp_new.xml inp.xml')
               INQUIRE (file='qfix', exist=l_exist)
               IF (l_exist) THEN
                  OPEN (2, file='qfix')
                  WRITE (2, *) "F"
                  CLOSE (2)
                  PRINT *, "qfix set to F"
               ENDIF
               CALL mixing_history_reset(fmpi)
            ENDIF
            call mixing_history_reset(fmpi)
         ENDIF
         CALL juDFT_end(" GEO new inp.xml created ! ", fmpi%irank)
      END SUBROUTINE priv_geo_end

   END SUBROUTINE fleur_execute
END MODULE m_fleur
