!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfprep_bulk
   !Reading of a FLEUR setup (bulk supercell or a written layer) and the
   !generation of the bulk potentials from the converged density.
   USE m_juDFT
   USE m_constants
   USE m_types_fleurinput
   USE m_types_wannier
   USE m_types_xcpot
   USE m_types_xcpot_inbuild
   USE m_types_sphhar
   USE m_types_stars
   USE m_types_mpi
   USE m_types_nococonv
   USE m_types_potden
   USE m_types_misc
   USE m_types_hybdat
   USE m_types_forcetheo
   USE m_types_lapw, ONLY: lapw_dim_nbasfcn
   USE m_fleurinput_read_xml
   USE m_fleurinput_postprocess
   USE m_make_xcpot
   USE m_make_sphhar
   USE m_make_stars
   USE m_make_forcetheo
   USE m_lapwdim
   USE m_convn
   USE m_cdn_io
   USE m_vgen
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_gfprep_struct, gfprep_read_struct, gfprep_bulk_init

   !a FLEUR setup as needed by gfprep
   TYPE t_gfprep_struct
      TYPE(t_fleurinput)  :: fi
      TYPE(t_wann)        :: wann
      CLASS(t_xcpot), ALLOCATABLE :: xcpot
      TYPE(t_sphhar)      :: sphhar
      TYPE(t_stars)       :: stars
   END TYPE t_gfprep_struct

CONTAINS

   SUBROUTINE gfprep_read_struct(filename_add, outxmlFileID, s, l_stars)
      !Reads <filename_add>inp.xml and sets up the lattice harmonics and
      !(unless l_stars=F) the stars exactly as fleur/gfleur do. The step
      !function is not generated.
      CHARACTER(LEN=*), INTENT(IN)       :: filename_add
      INTEGER, INTENT(IN)                :: outxmlFileID
      TYPE(t_gfprep_struct), INTENT(OUT) :: s
      LOGICAL, INTENT(IN), OPTIONAL      :: l_stars

      TYPE(t_kpts), ALLOCATABLE :: kptsArray(:)
      CHARACTER(LEN=40)  :: kptsSelection(3)
      CHARACTER(LEN=200) :: filename_add_loc

      filename_add_loc = filename_add
      ALLOCATE (t_xcpot_inbuild::s%xcpot)
      CALL fleurinput_read_xml(outxmlFileID, filename_add_loc, cell=s%fi%cell, sym=s%fi%sym, &
                               atoms=s%fi%atoms, input=s%fi%input, noco=s%fi%noco, &
                               vacuum=s%fi%vacuum, field=s%fi%field, sliceplot=s%fi%sliceplot, &
                               banddos=s%fi%banddos, mpinp=s%fi%mpinp, hybinp=s%fi%hybinp, &
                               coreSpecInput=s%fi%coreSpecInput, wann=s%wann, xcpot=s%xcpot, &
                               forcetheo_data=s%fi%forcetheo_data, kpts=s%fi%kpts, &
                               kptsSelection=kptsSelection, kptsArray=kptsArray, &
                               enparaXML=s%fi%enparaXML, gfinp=s%fi%gfinp, hub1inp=s%fi%hub1inp, &
                               dfpt=s%fi%dfpt)
      CALL fleurinput_postprocess(s%fi%cell, s%fi%sym, s%fi%atoms, s%fi%input, s%fi%noco, &
                                  s%fi%vacuum, s%fi%banddos, s%fi%hybinp, s%xcpot, &
                                  s%fi%kpts, s%fi%gfinp)
      CALL make_sphhar(.FALSE., s%fi%atoms, s%sphhar, s%fi%sym, s%fi%cell)
      IF (PRESENT(l_stars)) THEN
         IF (.NOT. l_stars) RETURN
      END IF
      s%stars%gmax = s%fi%input%gmax
      IF (ABS(s%fi%input%gmaxz) >= 1.0E-8) s%stars%gmaxz = s%fi%input%gmaxz
      CALL s%stars%dim(s%fi%sym, s%fi%cell, s%fi%input%film)
      CALL s%stars%init(s%fi%cell, s%fi%sym, s%fi%input%film, s%fi%input%rkmax)
   END SUBROUTINE gfprep_read_struct

   SUBROUTINE gfprep_bulk_init(fmpi, bulkdir, outxmlFileID, b, den, vTot, vCoul)
      !Reads the bulk supercell and its converged density from bulkdir and
      !generates the total and the Coulomb potential.
      TYPE(t_mpi), INTENT(INOUT)           :: fmpi
      CHARACTER(LEN=*), INTENT(IN)         :: bulkdir
      INTEGER, INTENT(IN)                  :: outxmlFileID
      TYPE(t_gfprep_struct), INTENT(OUT)   :: b
      TYPE(t_potden), INTENT(OUT)          :: den, vTot, vCoul

      TYPE(t_nococonv) :: nococonv
      TYPE(t_results)  :: results
      TYPE(t_hybdat)   :: hybdat
      TYPE(t_potden)   :: vx, vxc, exc, EnergyDen
      CLASS(t_forcetheo), ALLOCATABLE :: forcetheo
      INTEGER :: nbasfcn, archiveType
      REAL    :: ef, dist
      LOGICAL :: l_qfix, l_exist

      CALL gfprep_read_struct(bulkdir, outxmlFileID, b, l_stars=.FALSE.)

      IF (b%fi%input%film) CALL juDFT_error("gfprep needs a bulk supercell", &
                                            hint="the input in "//TRIM(bulkdir)//" is a film", calledby="gfprep_bulk_init")
      IF (ANY(ABS(b%fi%cell%amat(1:2, 3)) > 1.0E-8) .OR. ANY(ABS(b%fi%cell%amat(3, 1:2)) > 1.0E-8)) then
        CALL juDFT_error("The third lattice vector must be perpendicular to the first two", &
                          hint="gfprep stacks the layers along z", calledby="gfprep_bulk_init")
      end if
      IF (b%fi%noco%l_noco) CALL juDFT_error("Non-collinear bulk calculations are not supported by gfprep", &
                                            calledby="gfprep_bulk_init")

      !remaining setup as in fleur_init (the step function is needed by vgen)
      CALL make_xcpot(fmpi, b%xcpot, b%fi%atoms, b%fi%input)
      CALL nococonv%init(b%fi%noco)
      CALL nococonv%init_ss(b%fi%noco, b%fi%atoms)
      CALL make_stars(b%stars, b%fi%sym, b%fi%atoms, b%fi%vacuum, b%sphhar, b%fi%input, &
                      b%fi%cell, b%fi%noco, fmpi)
      CALL make_forcetheo(b%fi%forcetheo_data, b%fi%cell, b%fi%sym, b%fi%atoms, forcetheo)
      CALL b%fi%dfpt%init(b%fi%cell, b%fi%input)
      CALL lapw_dim(b%fi%kpts, b%fi%cell, b%fi%input, b%fi%noco, nococonv, forcetheo, &
                    b%fi%atoms, nbasfcn, b%fi%dfpt)
      CALL b%fi%input%init(b%fi%noco, b%fi%hybinp%l_hybrid, b%fi%sym%invs, &
                           b%fi%atoms%n_denmat, b%fi%atoms%n_hia, lapw_dim_nbasfcn)
      CALL b%fi%kpts%init(b%fi%sym, b%fi%input%film, .FALSE., .FALSE.)
      CALL convn(.FALSE., b%fi%atoms, b%stars)
      CALL results%init(b%fi%input, b%fi%atoms, b%fi%kpts, b%fi%noco)

      !converged density of the bulk
      INQUIRE (FILE=TRIM(bulkdir)//"cdn.hdf", EXIST=l_exist)
      IF (.NOT. l_exist) CALL juDFT_error("No density found: "//TRIM(bulkdir)//"cdn.hdf", &
                                          calledby="gfprep_bulk_init")
      CALL den%init(b%stars, b%fi%atoms, b%sphhar, b%fi%vacuum, b%fi%noco, b%fi%input%jspins, POTDEN_TYPE_DEN)
      archiveType = CDN_ARCHIVE_TYPE_CDN1_const
      CALL readDensity(b%stars, b%fi%noco, b%fi%vacuum, b%fi%atoms, b%fi%cell, b%sphhar, b%fi%input, &
                       b%fi%sym, archiveType, CDN_INPUT_DEN_const, 0, ef, dist, l_qfix, den, &
                       inFilename=TRIM(bulkdir)//"cdn")
      WRITE (oUnit, "(a,es12.4)") "gfprep: bulk density read, last distance:", dist

      !potentials
      CALL vTot%init(b%stars, b%fi%atoms, b%sphhar, b%fi%vacuum, b%fi%noco, b%fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL vCoul%init(b%stars, b%fi%atoms, b%sphhar, b%fi%vacuum, b%fi%noco, b%fi%input%jspins, POTDEN_TYPE_POTCOUL)
      CALL vx%init(b%stars, b%fi%atoms, b%sphhar, b%fi%vacuum, b%fi%noco, b%fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL vxc%init(b%stars, b%fi%atoms, b%sphhar, b%fi%vacuum, b%fi%noco, b%fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL exc%init(b%stars, b%fi%atoms, b%sphhar, b%fi%vacuum, b%fi%noco, b%fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL vgen(hybdat, b%fi%field, b%fi%input, b%xcpot, b%fi%atoms, b%sphhar, b%stars, b%fi%vacuum, &
                b%fi%sym, b%fi%cell, b%fi%sliceplot, fmpi, results, b%fi%noco, nococonv, EnergyDen, &
                den, vTot, vx, vCoul, vxc, exc)
   END SUBROUTINE gfprep_bulk_init

END MODULE m_gfprep_bulk
