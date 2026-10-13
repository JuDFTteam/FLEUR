!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfprep
   !gfprep: converts a converged bulk supercell (stacked along z) into the
   !layered setup of the Green-function embedding code gfleur.
   !
   !Input:  <bulkdir>inp.xml, <bulkdir>cdn.hdf, optional gfprep_inp
   !Output: <prefix><i>_inp.xml, gf_inp (template), gf_setup.hdf,
   !        gf_pot.hdf, gf_cdn.hdf, gfprep_planes.txt, gfprep.out
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_juDFT
   USE m_constants
   USE m_types_mpi
   USE m_types_potden
   USE m_types_kpts
   USE m_hdf_tools, ONLY: hdf_init
   USE m_check_arguments, ONLY: new_argument
   USE m_xmlOutput, ONLY: startFleur_XMLOutput, getXMLOutputUnitNumber
   USE m_gfprep_inp
   USE m_gfprep_bulk
   USE m_gfprep_layer
   USE m_gfprep_potconv
   USE m_gfprep_setup
   USE m_gfprep_write
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gfprep_run

CONTAINS

   SUBROUTINE gfprep_run()
      TYPE(t_mpi)            :: fmpi
      TYPE(t_gfprep_inp)     :: inp
      TYPE(t_gfprep_struct)  :: b, l
      TYPE(t_potden)         :: den, vTot, vCoul
      TYPE(t_layergeo), ALLOCATABLE     :: geo(:)
      TYPE(t_gfprep_layer), ALLOCATABLE :: lays(:)
      TYPE(t_kpts)           :: kpts
      COMPLEX, ALLOCATABLE   :: vpw(:, :), cpw(:), dpw(:, :)
      REAL, ALLOCATABLE      :: vmt(:, :, :, :), dmt(:, :, :, :)
      INTEGER :: i, js, jspins, n, ierr, outxml, kgrid(2)
      REAL    :: dv
      CHARACTER(LEN=40) :: name
      CHARACTER(LEN=100) :: xml_add

#ifdef CPP_MPI
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, fmpi%isize, ierr)
      IF (fmpi%isize > 1) CALL juDFT_error("gfprep must be run on a single MPI rank", calledby="gfprep")
      fmpi%mpi_comm = MPI_COMM_SELF
#else
      fmpi%mpi_comm = 1
#endif
      fmpi%irank = 0; fmpi%isize = 1
      ierr = 0
      CALL hdf_init()
      !queried by the inp.xml writer
      CALL new_argument(0, "-noco", "write out noco parameters into the layer inp.xml", "")
      CALL new_argument(0, "-greensf", "write out green's function parameters into the layer inp.xml", "")
      CALL new_argument(0, "-dfpt", "write out dfpt input into the layer inp.xml", "")
      OPEN (oUnit, FILE="gfprep.out", FORM="formatted", STATUS="replace")
      xml_add = "gfprep_"
      CALL startFleur_XMLOutput(xml_add)
      outxml = getXMLOutputUnitNumber()

      CALL gfprep_read_inp(inp)
      CALL priv_remove_outputs()

      !bulk: input, density and potentials
      CALL timestart("bulk")
      CALL gfprep_bulk_init(fmpi, TRIM(inp%bulkdir), outxml, b, den, vTot, vCoul)
      CALL timestop("bulk")
      jspins = b%fi%input%jspins

      !optional constant shift of the potential
      dv = inp%vshift
      IF (ABS(dv - 999.0) < 1.0E-6) dv = -SUM(REAL(vTot%pw(1, :jspins)))/jspins
      IF (dv /= 0.0) THEN
         WRITE (oUnit, "(a,f12.6)") "gfprep: potential shifted by ", dv
         WRITE (*, "(a,f12.6)") "gfprep: potential shifted by ", dv
         DO js = 1, jspins
            vTot%pw(1, js) = vTot%pw(1, js) + dv
            vCoul%pw(1, js) = vCoul%pw(1, js) + dv
            DO n = 1, b%fi%atoms%ntype
               ASSOCIATE (jri => b%fi%atoms%jri(n))
                  vTot%mt(:jri, 0, n, js) = vTot%mt(:jri, 0, n, js) + dv*b%fi%atoms%rmsh(:jri, n)
                  vCoul%mt(:jri, 0, n, js) = vCoul%mt(:jri, 0, n, js) + dv*sfp_const
               END ASSOCIATE
            END DO
         END DO
      END IF

      !layer geometry and structures
      CALL gfprep_geometry(inp, b%fi%atoms, b%fi%cell, geo)
      ALLOCATE (lays(SIZE(geo)))
      DO i = 1, SIZE(geo)
         CALL gfprep_make_layer(b, geo(i), inp%gmax_scale, lays(i))
      END DO
      kgrid = inp%kgrid
      IF (ANY(kgrid < 1)) kgrid = b%fi%kpts%nkpt3(1:2)
      CALL gfprep_kpts(lays, kgrid, .NOT. b%fi%noco%l_soc, kpts)

      !per layer: inp.xml, re-read as gfleur does, convert and write
      DO i = 1, SIZE(geo)
         WRITE (name, "(a,i0,a)") TRIM(inp%prefix), i, "_"
         CALL timestart("layer "//TRIM(name))
         CALL gfprep_write_inpxml(TRIM(name)//"inp.xml", b, lays(i), kpts)
         CALL gfprep_read_struct(TRIM(name), 0, l)

         ALLOCATE (vpw(l%stars%ng3, jspins), cpw(l%stars%ng3), dpw(l%stars%ng3, jspins))
         ALLOCATE (vmt(l%fi%atoms%jmtd, 0:l%sphhar%nlhd, l%fi%atoms%ntype, jspins))
         ALLOCATE (dmt(l%fi%atoms%jmtd, 0:l%sphhar%nlhd, l%fi%atoms%ntype, jspins))
         DO js = 1, jspins
            CALL gfprep_pwconv(b%stars, b%fi%cell, l%stars, l%fi%cell, geo(i)%zc, geo(i)%d, vTot%pw(:, js), vpw(:, js))
            CALL gfprep_pwconv(b%stars, b%fi%cell, l%stars, l%fi%cell, geo(i)%zc, geo(i)%d, den%pw(:, js), dpw(:, js))
         END DO
         CALL gfprep_pwconv(b%stars, b%fi%cell, l%stars, l%fi%cell, geo(i)%zc, geo(i)%d, vCoul%pw(:, 1), cpw)
         CALL gfprep_mtconv(b, l, lays(i)%typmap, lays(i)%atmap, vTot%mt(:, :, :, :jspins), vmt)
         CALL gfprep_mtconv(b, l, lays(i)%typmap, lays(i)%atmap, den%mt(:, :, :, :jspins), dmt)
         CALL gfprep_check("potential", b, l, i, geo(i)%zc, geo(i)%c, lays(i)%typmap, lays(i)%atmap, &
                           vTot%pw(:, 1), vpw(:, 1), vTot%mt, vmt)
         CALL gfprep_check("density", b, l, i, geo(i)%zc, geo(i)%c, lays(i)%typmap, lays(i)%atmap, &
                           den%pw(:, 1), dpw(:, 1), den%mt, dmt)

         IF (inp%spinflip(i) .AND. jspins == 2) THEN
            WRITE (oUnit, "(a,i0)") "gfprep: spins flipped in layer ", i
            vpw = vpw(:, [2, 1]); dpw = dpw(:, [2, 1])
            vmt = vmt(:, :, :, [2, 1]); dmt = dmt(:, :, :, [2, 1])
         END IF
         CALL gfprep_write_potden(i, l, vmt, vpw, cpw, dmt, dpw)
         CALL gfprep_write_setup(i, geo(i), b%fi%atoms, b%fi%cell)
         DEALLOCATE (vpw, cpw, dpw, vmt, dmt)
         CALL timestop("layer "//TRIM(name))
      END DO

      CALL gfprep_write_gfinp(inp%prefix, geo)
      WRITE (*, "(a,i0,a)") "gfprep: ", SIZE(geo), " layers written, see gfprep_planes.txt and gf_inp"
   END SUBROUTINE gfprep_run

   SUBROUTINE priv_remove_outputs()
      !Results of an earlier run would be mixed with the new ones
      CHARACTER(LEN=16), PARAMETER :: files(5) = [CHARACTER(LEN=16) :: "gf_setup.hdf", "gf_pot.hdf", &
                                                  "gf_cdn.hdf", "gf_steps.hdf", "gfprep_planes.txt"]
      INTEGER :: i
      LOGICAL :: l_exist
      DO i = 1, SIZE(files)
         INQUIRE (FILE=TRIM(files(i)), EXIST=l_exist)
         IF (.NOT. l_exist) CYCLE
         WRITE (*, *) "gfprep: removing old ", TRIM(files(i))
         OPEN (97, FILE=TRIM(files(i)))
         CLOSE (97, STATUS="delete")
      END DO
   END SUBROUTINE priv_remove_outputs

END MODULE m_gfprep

PROGRAM gfprep
   USE m_gfprep
   USE m_fleur_jobs, ONLY: fleur_job_init
   USE m_juDFT
   IMPLICIT NONE
   LOGICAL :: l_mpi_multithreaded
   CALL fleur_job_init(l_mpi_multithreaded)
   CALL gfprep_run()
   CALL juDFT_end("gfprep done")
END PROGRAM gfprep
