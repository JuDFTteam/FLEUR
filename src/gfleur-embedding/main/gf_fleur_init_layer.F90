!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_fleur_init_layer
   !Initializes one layer of the embedding calculation from its own
   !inp.xml file (<prefix><ilayer>_inp.xml) using the standard FLEUR
   !machinery - a per-layer "fleur_init light". Replaces the old
   !gf_setup.hdf handoff.
   !
   !Deliberately omitted compared to fleur_init: out.xml handling (one
   !stream for all layers, managed by gf_init), setupMPI (gfleur has
   !its own k x layer x energy decomposition), wannier/hybrid/
   !forcetheo modes, deleteDensities, usage data.
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_juDFT
   IMPLICIT NONE
CONTAINS
   SUBROUTINE gf_fleur_init_layer(fmpi, filename_add, outxmlFileID, ld)
      USE m_types
      USE m_constants
      USE m_gf_types
      USE m_fleurinput_read_xml
      USE m_fleurinput_postprocess
      USE m_fleurinput_mpi_bc
      USE m_make_xcpot
      USE m_make_sphhar
      USE m_make_stars
      USE m_make_forcetheo
      USE m_lapwdim
      USE m_convn
      USE m_cdn_io, ONLY: storeStructureIfNew
      USE m_gf_enpara, ONLY: gf_apply_enpara_atoms
      IMPLICIT NONE
      TYPE(t_mpi), INTENT(INOUT)     :: fmpi
      CHARACTER(len=*), INTENT(IN)   :: filename_add
      INTEGER, INTENT(IN)            :: outxmlFileID
      TYPE(t_gflayer), INTENT(INOUT) :: ld

      TYPE(t_enparaXML)      :: enparaXML
      TYPE(t_forcetheo_data) :: forcetheo_data
      TYPE(t_wann)           :: wann
      CLASS(t_forcetheo), ALLOCATABLE :: forcetheo
      TYPE(t_kpts), ALLOCATABLE       :: kptsArray(:)
      CHARACTER(LEN=40)      :: kptsSelection(3)
      CHARACTER(len=100)     :: filename_add_loc
      INTEGER                :: nbasfcn

      filename_add_loc = filename_add

      ALLOCATE (t_xcpot_inbuild::ld%xcpot)
      !Only PE==0 reads the input and does basic postprocessing
      IF (fmpi%irank == 0) THEN
         CALL fleurinput_read_xml(outxmlFileID, filename_add_loc, cell=ld%fi%cell, sym=ld%fi%sym, &
                                  atoms=ld%fi%atoms, input=ld%fi%input, noco=ld%fi%noco, &
                                  vacuum=ld%fi%vacuum, field=ld%fi%field, sliceplot=ld%fi%sliceplot, &
                                  banddos=ld%fi%banddos, mpinp=ld%fi%mpinp, hybinp=ld%fi%hybinp, &
                                  coreSpecInput=ld%fi%coreSpecInput, wann=wann, xcpot=ld%xcpot, &
                                  forcetheo_data=forcetheo_data, kpts=ld%fi%kpts, &
                                  kptsSelection=kptsSelection, kptsArray=kptsArray, &
                                  enparaXML=enparaXML, gfinp=ld%fi%gfinp, hub1inp=ld%fi%hub1inp, &
                                  dfpt=ld%fi%dfpt)
         CALL fleurinput_postprocess(ld%fi%cell, ld%fi%sym, ld%fi%atoms, ld%fi%input, ld%fi%noco, &
                                     ld%fi%vacuum, ld%fi%banddos, ld%fi%hybinp, ld%xcpot, &
                                     ld%fi%kpts, ld%fi%gfinp)
      END IF
      !Distribute the input to all PE
      CALL fleurinput_mpi_bc(ld%fi%cell, ld%fi%sym, ld%fi%atoms, ld%fi%input, ld%fi%noco, &
                             ld%fi%vacuum, ld%fi%field, ld%fi%sliceplot, ld%fi%banddos, &
                             ld%fi%mpinp, ld%fi%hybinp, ld%fi%coreSpecInput, wann, ld%xcpot, &
                             forcetheo_data, ld%fi%kpts, enparaXML, ld%fi%gfinp, ld%fi%hub1inp, &
                             fmpi%mpi_comm, ld%fi%dfpt)
      !Remaining init is done using all PE
      CALL make_xcpot(fmpi, ld%xcpot, ld%fi%atoms, ld%fi%input)
      CALL ld%nococonv%init(ld%fi%noco)
      CALL ld%nococonv%init_ss(ld%fi%noco, ld%fi%atoms)
      CALL ld%enpara%init_enpara(ld%fi%atoms, ld%fi%input%jspins, ld%fi%input%film, enparaXML)
      !apply the optional enpara_atoms override (quantum numbers)
      CALL gf_apply_enpara_atoms(ld%fi%atoms, ld%fi%input%jspins, ld%enpara)
      CALL make_sphhar(fmpi%irank == 0, ld%fi%atoms, ld%sphhar, ld%fi%sym, ld%fi%cell)
      !Store structure data (has to be performed before calling make_stars)
      CALL storeStructureIfNew(ld%fi%input, ld%stars, ld%fi%atoms, ld%fi%cell, ld%fi%vacuum, &
                               ld%fi%sym, fmpi, ld%sphhar, ld%fi%noco)
      CALL make_stars(ld%stars, ld%fi%sym, ld%fi%atoms, ld%fi%vacuum, ld%sphhar, ld%fi%input, &
                      ld%fi%cell, ld%fi%noco, fmpi)
      CALL make_forcetheo(forcetheo_data, ld%fi%cell, ld%fi%sym, ld%fi%atoms, forcetheo)
      CALL ld%fi%dfpt%init(ld%fi%cell, ld%fi%input)
      CALL lapw_dim(ld%fi%kpts, ld%fi%cell, ld%fi%input, ld%fi%noco, ld%nococonv, forcetheo, &
                    ld%fi%atoms, nbasfcn, ld%fi%dfpt)
      !snapshot the global LAPW dimensions for this layer - they must be
      !restored via lapw%init_dim before per-layer kernels run
      ld%nvd = lapw_dim_nvd
      ld%nv2d = lapw_dim_nv2d
      ld%nbasfcn = nbasfcn
      CALL ld%fi%input%init(ld%fi%noco, ld%fi%hybinp%l_hybrid, ld%fi%sym%invs, &
                            ld%fi%atoms%n_denmat, ld%fi%atoms%n_hia, lapw_dim_nbasfcn)
      !gfleur always inverts the complex (zS - H - Sigma): force complex
      !storage even if the layer has inversion symmetry
      ld%fi%input%l_real = .FALSE.
      CALL ld%fi%kpts%init(ld%fi%sym, ld%fi%input%film, .FALSE., .FALSE.)
      CALL convn(fmpi%irank == 0, ld%fi%atoms, ld%stars)

      !potential/density containers of the layer (filled during the SCF)
      CALL ld%vTot%init(ld%stars, ld%fi%atoms, ld%sphhar, ld%fi%vacuum, &
                        ld%fi%noco, ld%fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL ld%cdn_new%init(ld%stars, ld%fi%atoms, ld%sphhar, ld%fi%vacuum, &
                           ld%fi%noco, ld%fi%input%jspins, POTDEN_TYPE_DEN)
      ALLOCATE (ld%qmtl_new(0:MAXVAL(ld%fi%atoms%lmax), ld%fi%atoms%ntype))
      ld%qmtl_new = 0.0

      IF (.NOT. ld%fi%input%film) &
         CALL juDFT_error("gfleur layers must be film calculations", &
                          hint="set film geometry in "//TRIM(filename_add)//"inp.xml", &
                          calledby="gf_fleur_init_layer")
   END SUBROUTINE gf_fleur_init_layer
END MODULE m_gf_fleur_init_layer
