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
    TYPE(t_cell), INTENT(INOUT)  ::cell

    TYPE(t_vacuum), INTENT(INOUT)::vacuum
    TYPE(t_input), INTENT(INOUT) ::input
    TYPE(t_stars), INTENT(INOUT) ::stars
    TYPE(t_xcpot_inbuild_nf), INTENT(INOUT) ::xcpot
    TYPE(t_noco), INTENT(INOUT)  ::noco
    TYPE(t_mpinp), INTENT(INOUT) ::mpinp
    TYPE(t_hybinp), INTENT(INOUT)::hybinp

    INTEGER :: n
    REAL    :: min_dtild
    !
    !input
    !
    input%delgau = input%tkb
    IF (input%jspins==1.AND.(noco%l_ss.OR.noco%l_noco)) CALL judft_error("You cannot run a non-collinear calculation with a single spin, set jspins=2")
    IF (noco%l_ss) noco%l_noco = .TRUE.
    !check for magnetism
    DO n = 1, atoms%ntype
       IF (ANY(atoms%econf(n)%occupation(:, 1) .NE. atoms%econf(n)%occupation(:, 2))) THEN
          IF (input%jspins==1) CALL judft_error("You cannot run set different occupations for the two spins and use jspins=1")
          input%jspins = 2
       ENDIF
    ENDDO
    IF (input%jspins==0) input%jspins=1

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
       input%ctail = .FALSE.
    END IF

    IF (input%rkmax == 0.0) THEN
       input%rkmax = MAXVAL(atoms%lmax/atoms%rmt)
       input%rkmax = round_to_deci(input%rkmax, 1)
    ENDIF
    IF (input%rkmax > 4.5) THEN
       PRINT *, "WARNING, large default rkmax has been reduced. Check input"
       input%rkmax = 4.5
    ENDIF
    IF (noco%l_ss) input%ctail = .FALSE.
    input%zelec = DOT_PRODUCT(atoms%econf(:)%valence_electrons, atoms%neq(:))
    !
    ! stars
    !
    stars%gmax = MERGE(stars%gmax, round_to_deci(3.0*input%rkmax, 1), stars%gmax > 0)
    input%gmax = stars%gmax

    !
    !xcpot
    !
    xcpot%gmaxxc = MERGE(xcpot%gmaxxc, round_to_deci(3.0*input%rkmax,1), xcpot%gmaxxc > 0)
    xcpot%gmaxxc = MIN(input%gmax,xcpot%gmaxxc)

    xcpot%l_inbuild = .TRUE.
    IF (xcpot%icorr == 0) THEN
       IF (hybinp%l_hybrid) THEN
          xcpot%inbuild_name = "pbe0"
       ELSE
          xcpot%inbuild_name = "pbe"
       ENDIF
       CALL xcpot%init(atoms%ntype)
    ENDIF
    !
    !vacuum
    !
    IF (.NOT. input%film) THEN
       vacuum%dvac = cell%amat(3, 3)
    ELSE
       IF (vacuum%dvac<=0.0) THEN
          min_dtild=0.0
          DO n=1,atoms%ntype
             min_dtild=MAX(MAXVAL(ABS(atoms%pos(3,SUM(atoms%neq(:n-1))+1:SUM(atoms%neq(:n))))+atoms%rmt(n)),min_dtild)
          ENDDO
          IF (ABS(vacuum%dvac)<1E-4)THEN
             vacuum%dvac=2*min_dtild+0.2
          ELSE
             vacuum%dvac=2*min_dtild-vacuum%dvac
          ENDIF
          vacuum%dvac = round_to_deci(vacuum%dvac, 2)
          cell%amat(3,3)=vacuum%dvac+2.5
       ENDIF

    ENDIF
    vacuum%nvac = 2
    IF (sym%zrfs .OR. sym%invs) vacuum%nvac = 1
    !IF (oneD%odd%d1) vacuum%nvac = 1

    !
    !noco
    !
    ALLOCATE (noco%l_relax(atoms%ntype))
    ALLOCATE ( noco%alph_inp(atoms%ntype), noco%beta_inp(atoms%ntype))
    noco%qss_inp = MERGE(noco%qss_inp, [0.0, 0.0, 0.0], noco%l_ss)
    noco%l_relax(:) = .FALSE.
    noco%alph_inp(:) = 0.0
    noco%beta_inp(:) = 0.0

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
    CALL input%init(noco,hybinp%l_hybrid)
    hybinp%bands1 = CEILING(0.75*input%gw_neigd)
  END SUBROUTINE make_defaults
END MODULE m_make_defaults
