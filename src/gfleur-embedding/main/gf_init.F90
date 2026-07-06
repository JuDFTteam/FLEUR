!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_init
   !Setup of the embedding calculation. Replaces the old gf_setup:
   !instead of reading a gf_setup.hdf dumped by a modified FLEUR, each
   !layer is initialized from its own inp.xml via gf_fleur_init_layer;
   !gf_inp only carries the GF-specific control parameters.
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gf_init
CONTAINS
   SUBROUTINE gf_init(gmpi, embinp, gfmix, layers, ld, gkpts)
      USE m_types
      USE m_constants
      USE m_gf_types
      USE m_gf_out
      USE m_gf_version
      USE m_gf_iogfinp, ONLY: gf_rgfinp
      USE m_gf_fleur_init_layer
      USE m_gf_energies, gf_energies_init => init
      USE m_gf_mpi_groups
      USE m_gf_writetrans
      USE m_gaunt, ONLY: gaunt_init
      USE m_xmlOutput
#ifdef CPP_HDF
      USE m_hdf_tools, ONLY: hdf_init
#endif
      IMPLICIT NONE
      TYPE(t_gfmpi), INTENT(INOUT)  :: gmpi !fmpi%mpi_comm must be set
      TYPE(t_embinp), INTENT(OUT)   :: embinp
      TYPE(t_gfmix), INTENT(OUT)    :: gfmix
      TYPE(t_layers), INTENT(OUT)   :: layers
      TYPE(t_gflayer), ALLOCATABLE, INTENT(OUT) :: ld(:)
      TYPE(t_kpts), INTENT(OUT)     :: gkpts !the common 2D BZ k-points

      INTEGER :: i, outxmlFileID, ierr
      CHARACTER(len=100) :: filename_add, xml_add

#ifdef CPP_MPI
      CALL MPI_COMM_RANK(gmpi%fmpi%mpi_comm, gmpi%fmpi%irank, ierr)
      CALL MPI_COMM_SIZE(gmpi%fmpi%mpi_comm, gmpi%fmpi%isize, ierr)
#else
      gmpi%fmpi%irank = 0; gmpi%fmpi%isize = 1; gmpi%fmpi%mpi_comm = 1
#endif
      gmpi%pe0 = gmpi%fmpi%irank == 0
#ifdef CPP_HDF
      CALL hdf_init()
#endif
      CALL timestart("Setup")

      !open the out/inf files and say hello
      CALL gf_hello()
      CALL gf_out_newheader('SETUP of the Calculation')

      !one xml-output stream for the whole gfleur run; the per-layer
      !input reading documents itself into it
      outxmlFileID = -1
      IF (gmpi%pe0) THEN
         xml_add = ""
         CALL startFleur_XMLOutput(xml_add)
         outxmlFileID = getXMLOutputUnitNumber()
      END IF

      !Read the gf_inp-file (control parameters + layer decomposition)
      CALL timestart("reading setup")
      CALL gf_rgfinp(embinp, gfmix, layers)

      !Initialize each layer from its own inp.xml
      ALLOCATE (ld(layers%num_layers))
      DO i = 1, layers%num_layers
         WRITE (filename_add, "(a,i0,a)") TRIM(layers%prefix), i, "_"
         CALL timestart("init layer "//TRIM(filename_add))
         CALL gf_fleur_init_layer(gmpi%fmpi, TRIM(filename_add), outxmlFileID, ld(i))
         !distribute the GF-specific star cutoffs
         ld(i)%stars_gf%gmax_pot = embinp%gmax_pot
         ld(i)%stars_gf%gmax_decouple = embinp%gmax_decouple
         ld(i)%stars_gf%gmax_inp = ld(i)%stars%gmax
         CALL timestop("init layer "//TRIM(filename_add))
      END DO

      !gaunt tables sized for the largest lmax over all layers
      !(gaunt_init only acts on its first call)
      CALL gaunt_init(MAXVAL(ld(:)%fi%atoms%lmaxd) + 1)

      !cross-layer consistency checks
      CALL priv_check_layers(ld, layers)

      embinp%l_nogno = (.NOT. embinp%l_fullgreen .AND. (embinp%l_dos .OR. embinp%l_charge))

      !the 2D BZ is common to all layers - use the set of the first one
      gkpts = ld(1)%fi%kpts
      CALL timestop("reading setup")

      !read the energy contour
      CALL gf_out_newheader('Energy Grid for GFleur')
      CALL gf_energies_init(gmpi%pe0)

      !Setup the k x layer x energy parallelization
      CALL gf_setup_mpi_groups(gmpi, layers, gkpts, gf_noen())

      !sanity checks on the input switches
      CALL priv_setup_test(embinp, gfmix)

      IF (embinp%l_hdfio) THEN
         !open file for transmission data if needed
         IF (embinp%curr /= 0) CALL gf_writetrans_hdfopen(gmpi, gkpts, ld(1)%fi%sym, &
                                                          layers%num_layers, ld(1)%fi%input%jspins)
      END IF

      CALL gf_out_newheader('GF-SETUP done')
      CALL timestop("Setup")
   END SUBROUTINE gf_init

   !<-- S: priv_check_layers(ld,layers)
   SUBROUTINE priv_check_layers(ld, layers)
      !cross-layer consistency: all layers must share the in-plane cell,
      !the k-point set and the number of spins
      USE m_gf_types
      IMPLICIT NONE
      TYPE(t_gflayer), INTENT(IN) :: ld(:)
      TYPE(t_layers), INTENT(IN)  :: layers

      INTEGER :: i
      REAL, PARAMETER :: eps = 1.0E-8

      DO i = 2, layers%num_layers
         IF (ld(i)%fi%input%jspins /= ld(1)%fi%input%jspins) &
            CALL juDFT_error("jspins differs between the layers", calledby="gf_init")
         IF (ld(i)%fi%noco%l_noco .NEQV. ld(1)%fi%noco%l_noco) &
            CALL juDFT_error("l_noco differs between the layers", calledby="gf_init")
         IF (ANY(ABS(ld(i)%fi%cell%amat(1:2, 1:2) - ld(1)%fi%cell%amat(1:2, 1:2)) > eps)) &
            CALL juDFT_error("The in-plane lattices of the layers differ", calledby="gf_init")
         IF (ld(i)%fi%kpts%nkpt /= ld(1)%fi%kpts%nkpt) &
            CALL juDFT_error("The k-point sets of the layers differ", calledby="gf_init")
         IF (ANY(ABS(ld(i)%fi%kpts%bk(1:2, :) - ld(1)%fi%kpts%bk(1:2, :)) > eps)) &
            CALL juDFT_error("The k-point sets of the layers differ", calledby="gf_init")
      END DO
   END SUBROUTINE priv_check_layers
   !>

   !<-- S: priv_setup_test(embinp,gfmix)
   SUBROUTINE priv_setup_test(embinp, gfmix)
      !-----------------------------------------------
      !   Test the input for correct combination of switches
      !   More tests should be added here!!!
      !           (last modified: 2004-00-00) D. Wortmann
      !-----------------------------------------------
      USE m_constants, ONLY: oUnit
      USE m_gf_types
      IMPLICIT NONE
      TYPE(t_embinp), INTENT(IN)   :: embinp
      TYPE(t_gfmix), INTENT(INOUT) :: gfmix

      LOGICAL :: ok
      ok = .TRUE.

      !<-- Things that only should be warned about
      IF (.NOT. embinp%l_charge .AND. gfmix%iter > 0) THEN
         WRITE (oUnit, *) "WARNING, using iter>0 in non-selfconsistent calculation"
      END IF
      !>
      !<-- Tests that should disappear sometime :-)
      IF (embinp%npw /= 0) &
         CALL juDFT_error("Region II can only be treated analytically", calledby="gf_init")
      !>
      !<-- Test for CBS-mode
      IF (embinp%l_CBS) THEN
         gfmix%iter = 1
         IF (.NOT. embinp%l_gmat .AND. .NOT. embinp%l_tmat) ok = error( &
                                       "Check the tmat,gmat-switches for CBS-mode")
         IF (embinp%l_charge) ok = error("charge-generation cannot be used in CBS-MODE")
         IF (embinp%l_addemb) ok = error("l_addemb cannot be used in CBS-mode")
      END IF
      !>
      !<-- Test for sc
      IF (embinp%l_charge) THEN
         IF (.NOT. embinp%l_addemb) ok = error( &
                                    "Addemb must be specified for generation of charge")
      END IF
      !>
      !Test for current
      IF (embinp%curr /= 0) THEN
         IF (.NOT. embinp%l_addemb) ok = error( &
                                    "Addemb must be specified for calculation of conductances")
      END IF
      !<-- unsorted misc tests
      IF (embinp%l_IEC .AND. embinp%l_addemb) ok = error( &
                                              "in IEC mode you must set addemb = F!")
      !>

      IF (.NOT. ok) CALL juDFT_error("GF-setup:setup_test")

   CONTAINS
      !Error Message output!
      FUNCTION error(string)
         CHARACTER(LEN=*), INTENT(IN) :: string
         LOGICAL :: error
         WRITE (oUnit, *) "Setup Error in gf_init detected"
         WRITE (oUnit, *) string
         WRITE (*, *) "Setup Error in gf_init detected"
         WRITE (*, *) string
         error = .FALSE.
      END FUNCTION error
   END SUBROUTINE priv_setup_test
   !>
END MODULE m_gf_init
