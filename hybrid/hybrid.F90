!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_calc_hybinp
   USE m_judft

CONTAINS

   SUBROUTINE calc_hybinp(eig_id, mpdata, hybinp, kpts, atoms, input,  mpi, noco, cell, oneD, &
                          enpara, results, sym, xcpot, v, iterHF)

      USE m_types_hybdat
      USE m_types
      USE m_mixedbasis
      USE m_coulombmatrix
      USE m_hf_init
      USE m_hf_setup
      USE m_hsfock
      USE m_eig66_io
      USE m_io_hybinp

      IMPLICIT NONE

      TYPE(t_xcpot_inbuild), INTENT(IN)    :: xcpot
      TYPE(t_mpi), INTENT(IN)    :: mpi
      TYPE(t_oneD), INTENT(IN)    :: oneD
      type(t_mpdata), intent(inout) :: mpdata
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_potden), INTENT(IN)    :: v

      INTEGER, INTENT(INOUT) :: iterHF
      INTEGER, INTENT(IN)    :: eig_id

      ! local variables
      INTEGER           :: jsp, nk, err
      TYPE(t_hybdat)    :: hybdat
      type(t_lapw)      :: lapw
      LOGICAL           :: init_vex = .TRUE. !In first call we have to init v_nonlocal
      LOGICAL           :: l_zref
      character(len=999):: msg
      REAL, ALLOCATABLE :: eig_irr(:, :)

      ! open(7465, file="iter_translator.txt", position="append")
      ! write (7465,*) iter, iterHF
      ! close(7465)

      CALL timestart("hybinp code")
      INQUIRE (file="v_x.mat", exist=hybinp%l_addhf)
      CALL open_hybinp_io1( sym%invs)

      IF (kpts%nkptf == 0) THEN
         CALL judft_error("kpoint-set of full BZ not available", &
                          hint="to generate kpts in the full BZ you should specify a k-mesh in inp.xml")
      END IF

      !Check if new non-local potential shall be generated
      hybinp%l_subvxc = hybinp%l_hybrid .AND. (.NOT. xcpot%is_name("exx"))
      !If this is the first iteration loop we can not calculate a new non-local potential
      hybinp%l_calhf = (results%last_distance >= 0.0) .AND. (results%last_distance < input%minDistance)
      IF (.NOT. hybinp%l_calhf) THEN
         hybinp%l_subvxc = hybinp%l_subvxc .AND. hybinp%l_addhf
         CALL timestop("hybinp code")
         RETURN
      ENDIF

      results%te_hfex%core = 0

      !Check if we are converged well enough to calculate a new potential
      CALL open_hybinp_io1b( sym%invs)
      hybinp%l_addhf = .TRUE.

      !In first iteration allocate some memory
      IF (init_vex) THEN
         if(allocated(hybdat%ne_eig)) deallocate(hybdat%ne_eig)
         allocate(hybdat%ne_eig(kpts%nkpt), source=0)

         if(allocated(hybdat%nbands)) then
            deallocate(hybdat%nbands, stat=err, errmsg=msg)
            if(err /= 0) THEN
               write (*,*) "errorcode", err
               write (*,*) "errormessage", msg
            endif
         endif

         allocate(hybdat%nbands(kpts%nkpt), source=0)

         if(allocated(hybdat%nobd)) deallocate(hybdat%nobd)
         allocate(hybdat%nobd(kpts%nkptf, input%jspins), source=0)

         if(allocated(hybdat%nbasm)) deallocate(hybdat%nbasm)
         allocate(hybdat%nbasm(kpts%nkptf), source=0)

         if(allocated(hybdat%div_vv)) deallocate(hybdat%div_vv)
         allocate(hybdat%div_vv(input%neig, kpts%nkpt, input%jspins), source=0.0)
         init_vex = .FALSE.
      END IF

      hybinp%l_subvxc = (hybinp%l_subvxc .AND. hybinp%l_addhf)
      IF (.NOT. ALLOCATED(results%w_iks)) allocate(results%w_iks(input%neig, kpts%nkpt, input%jspins))

      IF (hybinp%l_calhf) THEN
         iterHF = iterHF + 1

         !Delete broyd files
         CALL system("rm -f broyd*")

         !check if z-reflection trick can be used

         l_zref = (sym%zrfs .AND. (SUM(ABS(kpts%bk(3, :kpts%nkpt))) < 1e-9) .AND. .NOT. noco%l_noco)

         CALL timestart("Preparation for hybinp functionals")
         !    CALL juDFT_WARN ("hybinp functionals not working in this version")

         !construct the mixed-basis
         CALL timestart("generation of mixed basis")
         write (*,*) "iterHF = ", iterHF
         CALL mixedbasis(atoms, kpts,  input, cell, xcpot, mpdata, hybinp, enpara, mpi, v, iterHF)
         CALL timestop("generation of mixed basis")

         CALL open_hybinp_io2(mpdata, hybinp, input, atoms, sym%invs)

         CALL coulombmatrix(mpi, atoms, kpts, cell, sym, mpdata, hybinp, xcpot)

         CALL hf_init(mpdata, hybinp, atoms, input,  hybdat)
         CALL timestop("Preparation for hybinp functionals")
         CALL timestart("Calculation of non-local HF potential")
         DO jsp = 1, input%jspins
            call timestart("HF_setup")
            CALL HF_setup(mpdata,hybinp, input, sym, kpts,  atoms, &
                          mpi, noco, cell, oneD, results, jsp, enpara, eig_id, &
                          hybdat, sym%invs, v%mt(:, 0, :, :), eig_irr)
            call timestop("HF_setup")

            DO nk = 1, kpts%nkpt
               !DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride
               CALL lapw%init(input, noco, kpts, atoms, sym, nk, cell, l_zref)
               CALL hsfock(nk, atoms, mpdata, hybinp, lapw,  kpts, jsp, input, hybdat, eig_irr, sym, cell, &
                           noco, results, MAXVAL(hybdat%nobd(:,jsp)), xcpot, mpi)
            END DO
         END DO
         CALL timestop("Calculation of non-local HF potential")
         CALL close_eig(eig_id)

      ENDIF
      CALL timestop("hybinp code")
   END SUBROUTINE calc_hybinp
END MODULE m_calc_hybinp
