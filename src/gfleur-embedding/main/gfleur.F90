!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfleur_driver
   !**********************************************************************
   !     Driver of the embedding Green-function code on the basis of
   !     the FLEUR code.
   !
   !     Daniel Wortmann 2000-2004
   !     Layer decomposition: Frank Freimuth, November 2007
   !
   !     CALL tree:
   !     gfleur
   !     |-gf_init             (per-layer init from inp.xml + gf_inp)
   !     |-gfleur_execute      (SCF / transport / CBS ... loops)
   !**********************************************************************
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gfleur_run
CONTAINS
   SUBROUTINE gfleur_run(l_mpi_multithreaded)
      USE m_types
      USE m_constants
      USE m_gf_types
      USE m_gf_init
      LOGICAL, INTENT(IN) :: l_mpi_multithreaded

      TYPE(t_gfmpi)   :: gmpi
      TYPE(t_embinp)  :: embinp
      TYPE(t_gfmix)   :: gfmix
      TYPE(t_layers)  :: layers
      TYPE(t_gflayer), ALLOCATABLE :: ld(:)
      TYPE(t_kpts)    :: gkpts

      gmpi%fmpi%l_mpi_multithreaded = l_mpi_multithreaded
#ifdef CPP_MPI
      gmpi%fmpi%mpi_comm = MPI_COMM_WORLD
#else
      gmpi%fmpi%mpi_comm = 1
#endif

      CALL gf_init(gmpi, embinp, gfmix, layers, ld, gkpts)

      CALL gfleur_execute(gmpi, embinp, gfmix, layers, ld, gkpts)
   END SUBROUTINE gfleur_run

   SUBROUTINE gfleur_execute(gmpi, embinp, gfmix, layers, ld, gkpts)
      !The self-consistency / transport / CBS loops of gfleur, transplanted
      !from the old PROGRAM gfleur (v25) onto the modern FLEUR kernels and
      !the per-layer t_gflayer containers.
      !
      !Structure (faithful to the original):
      !  iteration loop it
      !  |- pseudo-charge + layered interstitial Coulomb  (charge mode)
      !  |- gf_vgen: potential of every layer + enpara/tlmplm update
      !  |- spin loop jspin
      !  |  |- k loop nk
      !  |  |  |- layer loop: gf_apws + phasematrix + gf_gfleur_basic (H/S,T)
      !  |  |  |- gf_gfleur_propagate  (embedding self-energy propagation)
      !  |  |  |- layer loop: gf_gfleur_compose (G0 -> G, DOS, current, ...)
      !  |  |- gf_gencdn (charge mode)
      !  |- gf_mix / gf_totalmix (charge mode)
      !
      !ld(:) is the single source of truth; the few multi-layer routines
      !that still take plain arrays (Coulomb solver, mixing, DOS setup) get
      !temporary gathered arrays built from ld right before the call.
      USE m_types
      USE m_constants
      USE m_gf_types
      USE m_gf_optional, ONLY: gf_OPTIONAL
      USE m_gf_vgen, ONLY: gf_vgen
      USE m_gf_makepseudocharge, ONLY: gf_makepseudocharge
      USE m_gf_intcoul, ONLY: gf_makeintcoulpot
      USE m_gf_boundaryembpot, ONLY: gf_boundaryembpot
      USE m_gf_apws, ONLY: gf_apws
      USE m_gf_PhaseMatrix, ONLY: initPhaseMatrix
      USE m_gf_curvy2dprojector, ONLY: gf_curvy2dealloc
      USE m_gf_tlmplminit, ONLY: gf_tlmplminit
      USE m_gf_gfleur_basic, ONLY: gf_gfleur_basic
      USE m_gf_gfleur_propagate, ONLY: gf_gfleur_propagate
      USE m_gf_gfleur_compose, ONLY: gf_gfleur_compose
      USE m_gf_gencdn, ONLY: gf_gencdn
      USE m_gf_charge, ONLY: gf_charge, gf_charge_sum
      USE m_gf_totalmix, ONLY: gf_totalmix
      USE m_gf_mix, ONLY: gf_mix
      USE m_gf_dos, ONLY: gf_dos_init, gf_writedos
      USE m_gf_bandbending, ONLY: gf_bandbending
      USE m_gf_embpot_postprocess, ONLY: gf_embpot_spectral
      USE m_gf_CBS, ONLY: outCBS_openfile, outCBS_closefile
      USE m_gf_energies, ONLY: gf_noen
      USE m_gf_io2dmat, ONLY: io2dmat_finalize => finalize
      USE m_gf_writetrans, ONLY: gf_writetrans_hdfclose
#ifdef CPP_MPI
      USE mpi
#endif
      IMPLICIT NONE
      TYPE(t_gfmpi), INTENT(INOUT)  :: gmpi
      TYPE(t_embinp), INTENT(INOUT) :: embinp
      TYPE(t_gfmix), INTENT(INOUT)  :: gfmix
      TYPE(t_layers), INTENT(IN)    :: layers
      TYPE(t_gflayer), INTENT(INOUT), TARGET :: ld(:)
      TYPE(t_kpts), INTENT(IN)      :: gkpts

      INTEGER :: it, jspin, nspins, nk, nk_loop, layer, layer_loop, i
      INTEGER :: nl, neigd, ierr
      INTEGER :: chargelayers
      LOGICAL :: l_exist_cdn
      COMPLEX :: pot_aux(2, 3)
      REAL    :: vac_pot
      REAL, ALLOCATABLE :: qtot_el(:), qtot_nuc(:)
      !per-(layer) modern LAPW list, refilled per k by gf_apws; the GF
      !bookkeeping is held separately to avoid INOUT aliasing when a layer
      !is passed together with its own lapw_gf
      TYPE(t_lapw), ALLOCATABLE    :: lapw(:)
      TYPE(t_lapw_gf), ALLOCATABLE :: lapw_gf(:)
      !gathered plain-array views for the multi-layer array routines
      TYPE(t_atoms), ALLOCATABLE   :: atoms(:)
      TYPE(t_stars), ALLOCATABLE   :: stars(:)
      TYPE(t_sphhar), ALLOCATABLE  :: sphhar(:)
      TYPE(t_cell), ALLOCATABLE    :: cell(:)
      TYPE(t_enpara), ALLOCATABLE  :: enpara(:)
      TYPE(t_potden), ALLOCATABLE  :: potential(:)

      nl = layers%num_layers

      IF (embinp%l_surface) CALL juDFT_error( &
         "surface/vacuum embedding is not yet ported in this build", &
         calledby="gfleur_execute")

      IF (gmpi%pe0) THEN
         WRITE (oUnit, *)
         WRITE (oUnit, "(a)") "GFLEUR setup summary:"
         WRITE (oUnit, "(a,i0)") "  layers : ", nl
         WRITE (oUnit, "(a,i0)") "  k-points (2D BZ): ", gkpts%nkpt
         WRITE (oUnit, "(a,i0)") "  spins  : ", ld(1)%fi%input%jspins
         DO i = 1, nl
            WRITE (oUnit, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "  layer ", i, &
               ": ntype=", ld(i)%fi%atoms%ntype, " nat=", ld(i)%fi%atoms%nat, &
               " ng3=", ld(i)%stars%ng3, " nvd=", ld(i)%nvd
         END DO
      END IF

      !per-layer LAPW containers (dimensions come from gf_init)
      ALLOCATE (lapw(nl), lapw_gf(nl))
      DO layer = 1, nl
         lapw_gf(layer) = ld(layer)%lapw_gf
      END DO

      CALL gf_OPTIONAL(ld(1)%fi%input%jspins, ld(1)%fi%atoms, ld(1)%fi%input, &
                       ld(1)%nococonv, ld(1)%fi%cell, embinp, lapw(1), &
                       lapw_gf(1), ld(1)%fi%noco, gkpts, ld(1)%fi%sym)

      IF (embinp%l_charge) ALLOCATE (qtot_el(0:nl), qtot_nuc(0:nl))

      !<-- Self-consistency / analysis loop
      CALL timestart("iteration loop")
      DO it = 1, gfmix%iter
         IF (gmpi%pe0) WRITE (*, *) "Starting iteration: ", it, " of ", MAX(1, gfmix%iter)

         INQUIRE (FILE="gf_cdn.hdf", EXIST=l_exist_cdn)
         IF (.NOT. l_exist_cdn .AND. embinp%l_charge) THEN
            IF (gmpi%pe0) WRITE (oUnit, *) "WARNING, no gf_cdn.hdf but l_charge=T"
            IF (gfmix%iter > 1) CALL juDFT_error("No self-consistency without gf_cdn.hdf")
         END IF

         !<-- pseudo charge + layered interstitial Coulomb potential
         IF ((embinp%l_charge .OR. embinp%l_CBS) .AND. l_exist_cdn) THEN
            CALL timestart("gf_makepseudocharge")
            DO layer_loop = 1, gmpi%kl_LayerperPE
               layer = gmpi%kl_layers(layer_loop)
               IF (gmpi%k_kpts(1) == 1) &
                  CALL gf_makepseudocharge(layer, ld(1)%fi%input%jspins, ld(layer), layers, gmpi)
            END DO
            CALL timestop("gf_makepseudocharge")
         END IF

         IF (embinp%l_charge .AND. l_exist_cdn) THEN
            CALL timestart("gf_makeintcoulpot")
            DO layer = 1, nl
               ld(layer)%vTot%pw = 0.0
               ld(layer)%vTot%mt = 0.0
            END DO
#ifdef CPP_MPI
            CALL MPI_BARRIER(gmpi%fmpi%mpi_comm, ierr)
#endif
            IF (gmpi%pe0) THEN
               CALL priv_gather(ld, layers, atoms=atoms, stars=stars, cell=cell, potential=potential)
               CALL gf_makeintcoulpot(ld(1)%fi%input%jspins, layers, stars, gmpi%fmpi, &
                                      embinp, potential, vac_pot, atoms, cell, ld(1)%fi%sym)
               DO layer = 1, nl
                  ld(layer)%vTot%pw = potential(layer)%pw
               END DO
            END IF
            CALL timestop("gf_makeintcoulpot")
#ifdef CPP_MPI
            CALL MPI_BARRIER(gmpi%fmpi%mpi_comm, ierr)
#endif
         END IF
         !>

         !<-- Potential generation + energy-parameter / local Hamiltonian update
         IF (embinp%l_gmat) THEN
            CALL timestart("gf_vgen")
            IF (gmpi%pe0) WRITE (oUnit, *) "Potential Generation"
            CALL gf_vgen(gmpi, embinp, gfmix, layers, ld, .TRUE., &
                         ld(1)%fi%input%jspins, pot_aux)
            IF (embinp%l_CBS .AND. l_exist_cdn .AND. gmpi%fmpi%isize == 1) THEN
               CALL priv_gather(ld, layers, atoms=atoms, stars=stars, cell=cell)
               CALL gf_boundaryembpot(layers, stars, cell, atoms, embinp)
            END IF
            CALL timestop("gf_vgen")

            !energy parameters + radial functions + local Hamiltonian
            DO layer_loop = 1, gmpi%kl_LayerperPE
               layer = gmpi%kl_layers(layer_loop)
               CALL ld(layer)%enpara%update(gmpi%fmpi, ld(layer)%fi%atoms, &
                                            ld(layer)%fi%vacuum, ld(layer)%fi%input, &
                                            ld(layer)%vTot, ld(layer)%hub1data)
               CALL gf_tlmplminit(ld(layer), gmpi%fmpi)
            END DO
         END IF
         !>

         nspins = ld(1)%fi%input%jspins
         IF (embinp%l_IEC) nspins = 1
         IF (ld(1)%fi%noco%l_noco) nspins = 1

         !<-- Spin loop
         DO jspin = 1, nspins
            IF (embinp%l_dos) THEN
               CALL priv_gather(ld, layers, atoms=atoms, enpara=enpara, potential=potential)
               CALL gf_dos_init(layers, embinp, atoms, ld(1)%fi%cell, ld(1)%fi%sym, &
                                lapw(1), lapw_gf(1), ld(1)%fi%input%jspins, gkpts%nkpt, &
                                potential, enpara, ld(1)%fi%noco%l_noco)
            END IF

            IF (embinp%l_CBS) &
               CALL outCBS_openfile(jspin, ld(1)%fi%input%jspins, gkpts, SUM(layers%c), &
                                    2*ld(1)%nv2d, embinp%l_hdfio)

            !<-- k loop
            DO nk_loop = 1, gmpi%k_kperPE
               nk = gmpi%k_kpts(nk_loop)
               DO layer_loop = 1, gmpi%kl_layerperPE
                  layer = gmpi%kl_layers(layer_loop)
                  IF (gmpi%pe0) CALL priv_layer_progress('G0 for layer:', layer, &
                                                         layer == 1, layer_loop == gmpi%kl_layerperPE)

                  IF (embinp%l_nohelpregion) CALL gf_curvy2dealloc()

                  CALL timestart("gf_apws")
                  !restore this layer's global LAPW dimensions
                  CALL lapw(layer)%init_dim(ld(layer)%nvd, ld(layer)%nv2d, ld(layer)%nbasfcn)
                  CALL gf_apws(ld(1)%fi%input%jspins, jspin, gkpts%bk(:, nk), gmpi%pe0, &
                               ld(layer)%fi%atoms, ld(layer)%fi%input, ld(layer)%fi%sym, &
                               ld(layer)%fi%noco, ld(layer)%nococonv, embinp, &
                               ld(layer)%fi%cell, lapw(layer), lapw_gf(layer), layer)
                  CALL initPhaseMatrix(jspin, lapw(layer), lapw_gf(layer), &
                                       ld(layer)%fi%cell, embinp, ld(layer)%fi%noco%l_noco)
                  CALL timestop("gf_apws")

                  IF (embinp%l_gmat) THEN
                     CALL timestart("gf_gfleur_basic")
                     CALL gf_gfleur_basic(jspin, layer, nk, embinp, layers, ld(layer), &
                                          gmpi, gkpts%bk, lapw(layer), lapw_gf(layer), pot_aux)
                     CALL timestop("gf_gfleur_basic")
                  END IF
               END DO !layer

#ifdef CPP_MPI
               CALL MPI_BARRIER(gmpi%fmpi%mpi_comm, ierr)
#endif
               !<-- propagate the embedding self-energy through the stack
               IF (embinp%l_gmat .AND. (((embinp%l_dos .OR. embinp%l_charge) .AND. nl > 1) &
                                        .OR. (embinp%curr > 3))) THEN
                  CALL timestart("gf_gfleur_propagate")
                  layer = gmpi%kl_layers(1)
                  CALL gf_gfleur_propagate(layers, gmpi, lapw(layer), lapw_gf(layer), &
                                           embinp, nk, jspin, ld, gkpts%bk)
                  CALL timestop("gf_gfleur_propagate")
               END IF
#ifdef CPP_MPI
               CALL MPI_BARRIER(gmpi%fmpi%mpi_comm, ierr)
#endif
               IF (embinp%curr > 15) CYCLE

               !<-- second layer loop: G0 -> G, DOS, current, ...
               chargelayers = 1
               IF (embinp%l_charge .OR. embinp%l_dos) chargelayers = nl
               DO layer_loop = 1, MIN(chargelayers, gmpi%kl_layerperPE)
                  IF (embinp%l_gmat .AND. ((embinp%l_dos .OR. embinp%l_charge) &
                                           .OR. gmpi%kl_layers(layer_loop) == 1)) THEN
                     layer = gmpi%kl_layers(layer_loop)
                     IF (.NOT. (embinp%l_charge .OR. embinp%l_doslayer(layer))) CYCLE
                     IF (gmpi%pe0) CALL priv_layer_progress('G for layer:', layer, &
                                                            layer == 1, layer_loop == gmpi%kl_layerperPE)
                     CALL timestart("gf_gfleur_compose")
                     CALL gf_gfleur_compose(layer, ld(layer)%fi%noco, embinp, layers, &
                                            nk, jspin, ld(layer)%fi%sym, ld(layer)%fi%cell, &
                                            gmpi, lapw(layer), lapw_gf(layer), gkpts, pot_aux, &
                                            ld(layer)%cdn_new, ld(layer)%qmtl_new, &
                                            ld(layer)%fi%atoms, ld(layer)%stars, &
                                            ld(layer)%sphhar, ld(layer)%vTot%mt, ld(layer)%enpara)
                     CALL timestop("gf_gfleur_compose")
                  END IF
               END DO

               IF (gmpi%pe0) WRITE (*, "(a,i5,a,i5,a,i1,a,i1,a,i2,a,i2,a)") &
                  "Finished: k:(", nk, "/", gkpts%nkpt, ")  Spin: (", jspin, "/", &
                  ld(1)%fi%input%jspins, ") Iter: (", it, "/", gfmix%iter, ")"
               IF (embinp%l_embspectral) &
                  CALL gf_embpot_spectral(jspin, nk, embinp, layers, lapw(1), lapw_gf(1), ld(1)%fi%cell, gkpts)
            END DO !nk
            !>

            !<-- construct the new charge density of the layers
            neigd = 0
            IF (embinp%l_charge .AND. embinp%l_gmat) THEN
               qtot_el = 0.0
               qtot_nuc = 0.0
               DO layer_loop = 1, gmpi%kl_layerperPE
                  layer = gmpi%kl_layers(layer_loop)
                  CALL timestart("gf_gencdn")
                  CALL gf_gencdn(layer, jspin, .FALSE., ld(1)%fi%input%jspins, embinp, &
                                 ld(layer)%fi%input, ld(layer)%fi%atoms, ld(layer)%fi%cell, &
                                 ld(layer)%fi%sym, gkpts, ld(layer)%stars, ld(layer)%sphhar, &
                                 gmpi, ld(layer)%enpara, ld(layer)%vTot%mt, ld(layer)%cdn_new%pw, &
                                 ld(layer)%cdn_new%mt, ld(layer)%qmtl_new, neigd, &
                                 ld(1)%fi%noco%l_noco, qtot_el(1:), qtot_nuc(1:))
                  CALL timestop("gf_gencdn")
                  ld(layer)%qmtl_new = 0.0
               END DO
            END IF
            !>
            IF (embinp%l_CBS) CALL outCBS_closefile(ld(1)%fi%noco%l_noco .OR. &
                                                    (jspin == ld(1)%fi%input%jspins))
         END DO !spin
         !>

         !<-- mix the charge density
         IF (embinp%l_charge .AND. l_exist_cdn) THEN
            CALL timestart("gf_mix")
            CALL gf_charge_sum(embinp%l_surface, qtot_nuc, qtot_el)
            IF (gmpi%pe0) WRITE (oUnit, "('Total charge neutrality:',f9.4,'+',f9.4,'=',f9.4)") &
               -qtot_el(0), qtot_nuc(0), -qtot_el(0) + qtot_nuc(0)
            DO layer_loop = 1, gmpi%kl_LayerperPE
               layer = gmpi%kl_layers(layer_loop)
               IF (gmpi%k_kpts(1) == 1) THEN
                  WHERE (ABS(qtot_el) < 1E-10) qtot_el = 10E-10
                  CALL gf_charge(ld(1)%fi%input%jspins, gmpi, ld(layer)%stars, &
                                 ld(layer)%fi%atoms, ld(layer)%fi%cell, gfmix, embinp, &
                                 ld(layer)%fi%noco, ld(layer)%sphhar, ld(layer)%fi%sym, &
                                 ld(layer)%cdn_new%pw, ld(layer)%cdn_new%mt, layer, &
                                 qtot_nuc/qtot_el)
               END IF
               ld(layer)%cdn_new%pw = 0.0
               ld(layer)%cdn_new%mt = 0.0
               ld(layer)%qmtl_new = 0.0
            END DO
            IF (.NOT. embinp%l_potmix) THEN
               IF (embinp%l_totalmix) THEN
                  CALL priv_gather(ld, layers, atoms=atoms, stars=stars, sphhar=sphhar)
                  CALL gf_totalmix(ld(1)%fi%input%jspins, layers, atoms, sphhar, stars, &
                                   gfmix, embinp, gmpi)
               ELSE
                  CALL gf_mix(gmpi, embinp, gfmix, layers, ld, ld(1)%fi%input%jspins)
               END IF
            END IF
            CALL timestop("gf_mix")
         END IF
         !>
#ifdef CPP_MPI
         CALL MPI_BARRIER(gmpi%fmpi%mpi_comm, ierr)
#endif
         IF (gmpi%pe0) WRITE (*, *) "Finished IT:", it
      END DO !iteration
      CALL timestop("iteration loop")

      IF (nl == 1 .AND. embinp%l_CBS) &
         CALL gf_bandbending(lapw(1), lapw_gf(1), ld(1)%fi%input%jspins, gkpts%nkpt)

      !<-- output
      CALL timestart("gf_writedos")
      IF (embinp%l_dos) THEN
         CALL priv_gather(ld, layers, atoms=atoms, cell=cell)
         IF (ld(1)%fi%noco%l_noco) THEN
            CALL gf_writedos(gmpi, layers, atoms, cell, embinp, 3, gkpts%nkpt)
         ELSE
            CALL gf_writedos(gmpi, layers, atoms, cell, embinp, ld(1)%fi%input%jspins, gkpts%nkpt)
         END IF
      END IF
      CALL timestop("gf_writedos")

      CALL io2dmat_finalize()
      CALL gf_writetrans_hdfclose()
      IF (gmpi%pe0) CALL writetimes()
   END SUBROUTINE gfleur_execute

   !<-- S: priv_gather(ld,layers,...)
   SUBROUTINE priv_gather(ld, layers, atoms, stars, sphhar, cell, enpara, potential)
      !Build plain-array copies of the requested per-layer components for
      !the multi-layer routines that still take arrays. Only the arrays
      !that are actually present are (re)allocated and filled.
      USE m_types
      USE m_gf_types
      IMPLICIT NONE
      TYPE(t_gflayer), INTENT(IN)  :: ld(:)
      TYPE(t_layers), INTENT(IN)   :: layers
      TYPE(t_atoms), ALLOCATABLE, INTENT(INOUT), OPTIONAL   :: atoms(:)
      TYPE(t_stars), ALLOCATABLE, INTENT(INOUT), OPTIONAL   :: stars(:)
      TYPE(t_sphhar), ALLOCATABLE, INTENT(INOUT), OPTIONAL  :: sphhar(:)
      TYPE(t_cell), ALLOCATABLE, INTENT(INOUT), OPTIONAL    :: cell(:)
      TYPE(t_enpara), ALLOCATABLE, INTENT(INOUT), OPTIONAL  :: enpara(:)
      TYPE(t_potden), ALLOCATABLE, INTENT(INOUT), OPTIONAL  :: potential(:)
      INTEGER :: i, nl
      nl = layers%num_layers
      IF (PRESENT(atoms)) THEN
         IF (ALLOCATED(atoms)) DEALLOCATE (atoms)
         ALLOCATE (atoms(nl)); DO i = 1, nl; atoms(i) = ld(i)%fi%atoms; END DO
      END IF
      IF (PRESENT(stars)) THEN
         IF (ALLOCATED(stars)) DEALLOCATE (stars)
         ALLOCATE (stars(nl)); DO i = 1, nl; stars(i) = ld(i)%stars; END DO
      END IF
      IF (PRESENT(sphhar)) THEN
         IF (ALLOCATED(sphhar)) DEALLOCATE (sphhar)
         ALLOCATE (sphhar(nl)); DO i = 1, nl; sphhar(i) = ld(i)%sphhar; END DO
      END IF
      IF (PRESENT(cell)) THEN
         IF (ALLOCATED(cell)) DEALLOCATE (cell)
         ALLOCATE (cell(nl)); DO i = 1, nl; cell(i) = ld(i)%fi%cell; END DO
      END IF
      IF (PRESENT(enpara)) THEN
         IF (ALLOCATED(enpara)) DEALLOCATE (enpara)
         ALLOCATE (enpara(nl)); DO i = 1, nl; enpara(i) = ld(i)%enpara; END DO
      END IF
      IF (PRESENT(potential)) THEN
         IF (ALLOCATED(potential)) DEALLOCATE (potential)
         ALLOCATE (potential(nl)); DO i = 1, nl; potential(i) = ld(i)%vTot; END DO
      END IF
   END SUBROUTINE priv_gather
   !>

   !<-- S: priv_layer_progress(text,layer,first,last)
   SUBROUTINE priv_layer_progress(text, layer, first, last)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: text
      INTEGER, INTENT(IN)          :: layer
      LOGICAL, INTENT(IN)          :: first, last
      IF (first) WRITE (*, "(A)", advance='no') text
      IF (last) THEN
         WRITE (*, "(i3)") layer
      ELSE
         WRITE (*, "(i3)", advance='no') layer
      END IF
   END SUBROUTINE priv_layer_progress
   !>
END MODULE m_gfleur_driver

PROGRAM gfleur
   USE m_gfleur_driver
   USE m_fleur_jobs, ONLY: fleur_job_init
   USE m_juDFT
   IMPLICIT NONE
   LOGICAL :: l_mpi_multithreaded
   CALL fleur_job_init(l_mpi_multithreaded)
   CALL gfleur_run(l_mpi_multithreaded)
   CALL juDFT_end("GFLEUR done")
END PROGRAM gfleur
