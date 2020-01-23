!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_defaults
   USE m_juDFT
   IMPLICIT NONE
!---------------------------------------------------------------------
!  Check muffin tin radii and determine a reasonable choice for MTRs.
!  Derive also other parameters for the input file, to provide some
!  help in the out-file.                                        gb`02
!---------------------------------------------------------------------
CONTAINS

   SUBROUTINE make_defaults(atoms, sym, cell, vacuum, input, stars, &
                            xcpot, noco, mpinp, hybinp)
      USE m_types_atoms
      USE m_types_cell
      USE m_types_sym
      USE m_types_vacuum
      USE m_types_xcpot_inbuild_nofunction
      USE m_types_input
      USE m_types_stars
      USE m_types_noco
      USE m_types_mpinp
      USE m_types_hybinp
      USE m_juDFT

      TYPE(t_atoms), INTENT(IN)    ::atoms
      TYPE(t_sym), INTENT(IN)      ::sym
      TYPE(t_cell), INTENT(IN)     ::cell

      TYPE(t_vacuum), INTENT(INOUT)::vacuum
      TYPE(t_input), INTENT(INOUT) ::input
      TYPE(t_stars), INTENT(INOUT) ::stars
      TYPE(t_xcpot_inbuild_nf), INTENT(INOUT) ::xcpot
      TYPE(t_noco), INTENT(INOUT)  ::noco
      TYPE(t_mpinp), INTENT(INOUT) ::mpinp
      TYPE(t_hybinp), INTENT(INOUT)::hybinp

      integer :: n
      !
      !input
      !
      input%delgau = input%tkb
      IF (noco%l_ss) noco%l_noco = .TRUE.
      IF (noco%l_noco) input%jspins = 2
      !check for magnetism
      DO n = 1, atoms%ntype
         IF (ANY(atoms%econf(n)%occupation(:, 1) .NE. atoms%econf(n)%occupation(:, 2))) input%jspins = 2
      ENDDO

      IF (ANY(atoms%nlo(:) .NE. 0)) THEN
         input%ellow = -1.8
      ELSE
         input%ellow = -0.8
      ENDIF
      IF (input%film) THEN
         input%elup = 0.5
      ELSE
         input%elup = 1.0
      ENDIF
      IF (hybinp%l_hybrid) THEN
         input%ellow = input%ellow - 2.0
         input%elup = input%elup + 10.0
         input%minDistance = 1.0e-5
         input%ctail = .false.
      END IF

      IF (input%rkmax == 0.0) input%rkmax = MAXVAL(atoms%lmax/atoms%rmt)
      IF (input%rkmax > 4.5) THEN
         PRINT *, "WARNING, large default rkmax has been reduced. Check input"
         input%rkmax = 4.5
      ENDIF
      input%rkmax = round_to_deci(input%rkmax, 1)
      IF (noco%l_ss) input%ctail = .FALSE.
      input%zelec = dot_product(atoms%econf(:)%valence_electrons, atoms%neq(:))
      !
      ! stars
      !
      stars%gmax = merge(stars%gmax, 3.0*input%rkmax, stars%gmax > 0)
      stars%gmax = round_to_deci(stars%gmax, 1)
      input%gmax = stars%gmax

      !
      !xcpot
      !
      xcpot%gmaxxc = merge(xcpot%gmaxxc, 3.0*input%rkmax, xcpot%gmaxxc > 0)
      xcpot%gmaxxc = min(input%gmax,xcpot%gmaxxc)
      xcpot%gmaxxc = round_to_deci(xcpot%gmaxxc, 1)
      xcpot%l_inbuild = .true.
      if (xcpot%icorr == 0) THEN
         if (hybinp%l_hybrid) then
            xcpot%inbuild_name = "pbe0"
         else
            xcpot%inbuild_name = "pbe"
         endif
         call xcpot%init(atoms%ntype)
      endif
      !
      !vacuum
      !
      IF (.not. input%film) THEN
         vacuum%dvac = cell%amat(3, 3)
      Else
         vacuum%dvac = round_to_deci(vacuum%dvac, 2)
      ENDIF
      vacuum%nvac = 2
      IF (sym%zrfs .OR. sym%invs) vacuum%nvac = 1
      !IF (oneD%odd%d1) vacuum%nvac = 1

      !
      !noco
      !
      ALLOCATE (noco%l_relax(atoms%ntype), noco%b_con(2, atoms%ntype))
      ALLOCATE (noco%alphInit(atoms%ntype), noco%alph(atoms%ntype), noco%beta(atoms%ntype))
      noco%qss = MERGE(noco%qss, [0.0, 0.0, 0.0], noco%l_ss)
      noco%l_relax(:) = .FALSE.
      noco%alphInit(:) = 0.0
      noco%alph(:) = 0.0
      noco%beta(:) = 0.0
      noco%b_con(:, :) = 0.0

      !
      !hybinp
      !
      ALLOCATE (hybinp%lcutwf(atoms%ntype))
      ALLOCATE (hybinp%lcutm1(atoms%ntype))
      ALLOCATE (hybinp%select1(4, atoms%ntype))
      hybinp%lcutwf = atoms%lmax - atoms%lmax/10
      hybinp%lcutm1 = 4
      hybinp%select1(1, :) = 4
      hybinp%select1(2, :) = 0
      hybinp%select1(3, :) = 4
      hybinp%select1(4, :) = 2
      mpinp%g_cutoff = round_to_deci(input%rkmax - 0.5, 1)
      mpinp%linear_dep_tol = 1e-4
      call input%init(noco,hybinp%l_hybrid)
      hybinp%bands1 = ceiling(0.75*input%gw_neigd)
   END SUBROUTINE make_defaults
END MODULE m_make_defaults
