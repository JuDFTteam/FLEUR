!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_calc_hybrid

   USE m_judft

CONTAINS

   SUBROUTINE calc_hybrid(eig_id, hybrid, kpts, atoms, input, DIMENSION, mpi, noco, cell, oneD, &
                          enpara, results, sym, xcpot, v, iter, iterHF)

      USE m_types
      USE m_mixedbasis
      USE m_coulombmatrix
      USE m_hf_init
      USE m_hf_setup
      USE m_hsfock
      USE m_eig66_io
      USE m_io_hybrid

      IMPLICIT NONE

      TYPE(t_xcpot_inbuild), INTENT(IN)    :: xcpot
      TYPE(t_mpi), INTENT(IN)    :: mpi
      TYPE(t_dimension), INTENT(IN)    :: DIMENSION
      TYPE(t_oneD), INTENT(IN)    :: oneD
      TYPE(t_hybrid), INTENT(INOUT) :: hybrid
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_potden), INTENT(IN)    :: v

      INTEGER, INTENT(IN)    :: iter
      INTEGER, INTENT(INOUT) :: iterHF
      INTEGER, INTENT(IN)    :: eig_id

      ! local variables
      INTEGER           :: jsp, nk, nred
      TYPE(t_hybdat)    :: hybdat
      type(t_lapw)      :: lapw
      LOGICAL           :: init_vex = .TRUE. !In first call we have to init v_nonlocal
      LOGICAL           :: l_restart = .FALSE.
      LOGICAL           :: l_zref

      REAL              :: bkpt(3)
      REAL, ALLOCATABLE :: eig_irr(:, :)

      CALL timestart("Hybrid code")
      INQUIRE (file="v_x.mat", exist=hybrid%l_addhf)
      CALL open_hybrid_io1(DIMENSION, sym%invs)

      IF (kpts%nkptf == 0) THEN
         CALL judft_error("kpoint-set of full BZ not available", &
                          hint="to generate kpts in the full BZ you should specify a k-mesh in inp.xml")
      END IF

      !Check if new non-local potential shall be generated
      hybrid%l_subvxc = hybrid%l_hybrid .AND. (.NOT. xcpot%is_name("exx"))
      !If this is the first iteration loop we can not calculate a new non-local potential
      hybrid%l_calhf = (results%last_distance >= 0.0) .AND. (results%last_distance < input%minDistance)
      IF (.NOT. hybrid%l_calhf) THEN
         hybrid%l_subvxc = hybrid%l_subvxc .AND. hybrid%l_addhf
         CALL timestop("Hybrid code")
         RETURN
      ENDIF

      results%te_hfex%core = 0

      !Check if we are converged well enough to calculate a new potential
      CALL open_hybrid_io1b(DIMENSION, sym%invs)
      hybrid%l_addhf = .TRUE.

      !In first iteration allocate some memory
      IF (init_vex) THEN
         ALLOCATE (hybrid%ne_eig(kpts%nkpt), hybrid%nbands(kpts%nkpt), hybrid%nobd(kpts%nkptf))
         ALLOCATE (hybrid%nbasm(kpts%nkptf))
         ALLOCATE (hybrid%div_vv(DIMENSION%neigd, kpts%nkpt, input%jspins))
         init_vex = .FALSE.
      END IF

      hybrid%l_subvxc = (hybrid%l_subvxc .AND. hybrid%l_addhf)
      IF (.NOT. ALLOCATED(results%w_iks)) ALLOCATE (results%w_iks(merge(dimension%neigd*2,dimension%neigd,noco%l_soc), kpts%nkpt, input%jspins))

      IF (hybrid%l_calhf) THEN
         iterHF = iterHF + 1

         !Delete broyd files
         CALL system("rm -f broyd*")

         !check if z-reflection trick can be used

         l_zref = (sym%zrfs .AND. (SUM(ABS(kpts%bk(3, :kpts%nkpt))) < 1e-9) .AND. .NOT. noco%l_noco)

         CALL timestart("Preparation for Hybrid functionals")
         !    CALL juDFT_WARN ("Hybrid functionals not working in this version")

         !construct the mixed-basis
         CALL timestart("generation of mixed basis")
         CALL mixedbasis(atoms, kpts, dimension, input, cell, sym, xcpot, hybrid, enpara, mpi, v, l_restart)
         CALL timestop("generation of mixed basis")

         CALL open_hybrid_io2(hybrid, DIMENSION, atoms, sym%invs)

         CALL coulombmatrix(mpi, atoms, kpts, cell, sym, hybrid, xcpot, l_restart)

         CALL hf_init(hybrid, kpts, atoms, input, DIMENSION, hybdat, sym%invs)
         CALL timestop("Preparation for Hybrid functionals")
         CALL timestart("Calculation of non-local HF potential")
         DO jsp = 1, input%jspins
            call timestart("HF_setup")
            CALL HF_setup(hybrid, input, sym, kpts, dimension, atoms, mpi, noco, cell, oneD, results, jsp, enpara, eig_id, &
                          hybdat, iterHF, sym%invs, v%mt(:, 0, :, :), eig_irr)
            call timestop("HF_setup")

            DO nk = 1, kpts%nkpt
               !DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride
               CALL lapw%init(input, noco, kpts, atoms, sym, nk, cell, l_zref)
               CALL hsfock(nk, atoms, hybrid, lapw, DIMENSION, kpts, jsp, input, hybdat, eig_irr, sym, cell, &
                           noco, results, iterHF, MAXVAL(hybrid%nobd), xcpot, mpi)
            END DO
         END DO
         CALL timestop("Calculation of non-local HF potential")
         CALL close_eig(eig_id)

      ENDIF
      CALL timestop("Hybrid code")
   END SUBROUTINE calc_hybrid
END MODULE m_calc_hybrid
