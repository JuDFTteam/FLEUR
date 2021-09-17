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
       xcpot, noco, banddos, mpinp, hybinp)
    USE m_types_atoms
    USE m_types_cell
    USE m_types_sym
    USE m_types_vacuum
    USE m_types_xcpot_inbuild_nofunction
    USE m_types_input
    USE m_types_stars
    USE m_types_noco
    USE m_types_banddos
    USE m_types_mpinp
    USE m_types_hybinp
    USE m_juDFT

    TYPE(t_atoms), INTENT(INOUT) ::atoms
    TYPE(t_sym), INTENT(IN)      ::sym
    TYPE(t_cell), INTENT(INOUT)  ::cell

    TYPE(t_vacuum), INTENT(INOUT)::vacuum
    TYPE(t_input), INTENT(INOUT) ::input
    TYPE(t_stars), INTENT(INOUT) ::stars
    TYPE(t_xcpot_inbuild_nf), INTENT(INOUT) ::xcpot
    TYPE(t_noco), INTENT(INOUT)  ::noco
    TYPE(t_banddos), INTENT(INOUT) :: banddos
    TYPE(t_mpinp), INTENT(INOUT) ::mpinp
    TYPE(t_hybinp), INTENT(INOUT)::hybinp

    INTEGER :: n
    REAL    :: min_dtild
    !
    !input
    !
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

    IF (hybinp%l_hybrid) THEN
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
       IF (cell%primCellZ.NE.0.0) THEN
          IF(cell%amat(3,3) - (cell%primCellZ * NINT(cell%amat(3,3) / cell%primCellZ)).GT.1.0e-5) THEN
             WRITE(*,'(a,f15.8)') 'Supercell z: ', cell%amat(3,3)
             WRITE(*,'(a,f15.8)') 'primitive cell z: ', cell%primCellZ
             WRITE(*,'(a,f15.8)') 'Factor between supercell z and primitive cell z: ', cell%amat(3,3) / cell%primCellZ
             CALL juDFT_warn("supercell and primitive cell z dimensions are inconsistent.", calledby="make_defaults")
          END IF
          banddos%s_cell_z = NINT(cell%amat(3,3) / cell%primCellZ)
       END IF
    ELSE
       IF (vacuum%dvac<=abs(cell%amat(3,3))) THEN
          min_dtild=0.0
          DO n=1,atoms%ntype
             min_dtild=MAX(MAXVAL(ABS(atoms%pos(3,SUM(atoms%neq(:n-1))+1:SUM(atoms%neq(:n))))+atoms%rmt(n)),min_dtild)
          ENDDO
          IF (ABS(vacuum%dvac)<=abs(cell%amat(3,3)))THEN
             vacuum%dvac=2*min_dtild+0.2
          ELSE
             vacuum%dvac=2*min_dtild-vacuum%dvac
          ENDIF
          vacuum%dvac = round_to_deci(vacuum%dvac, 2)
          atoms%taual(3,:)=atoms%taual(3,:)*cell%amat(3,3)/(vacuum%dvac+2.5)
       ENDIF
       cell%amat(3,3)=vacuum%dvac+2.5
       IF(cell%primCellZ.NE.0.0) THEN
          cell%amat(3,3) = cell%primCellZ * CEILING(cell%amat(3,3) / cell%primCellZ)
          banddos%s_cell_z = NINT(cell%amat(3,3) / cell%primCellZ)
       END IF
    ENDIF
    vacuum%nvac = 2
    IF (sym%zrfs .OR. sym%invs) vacuum%nvac = 1
    !IF (oneD%odd%d1) vacuum%nvac = 1

    !
    !noco
    !

    ALLOCATE ( noco%alph_inp(atoms%ntype), noco%beta_inp(atoms%ntype))
    noco%qss_inp = MERGE(noco%qss_inp, [0.0, 0.0, 0.0], noco%l_ss)
    noco%alph_inp(:) = 0.0
    noco%beta_inp(:) = 0.0
    ALLOCATE(noco%l_constrained(atoms%ntype))
    noco%l_constrained(:) = .FALSE.
    ALLOCATE(noco%l_unrestrictMT(atoms%ntype))
    noco%l_unrestrictMT(:) = .FALSE.
    ALLOCATE(noco%l_alignMT(atoms%ntype))
    noco%l_alignMT(:)=.false.
    ALLOCATE(noco%mix_RelaxWeightOffD(atoms%ntype))
    noco%mix_RelaxWeightOffD(:)=1.0
    noco%mag_mixing_scheme=0

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
