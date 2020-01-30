module m_calc_cmt

contains
   ! subroutine calc_cmt(atoms, cell, input, noco, hybinp, hybdat, kpts, sym, oneD,&
   !                     zmat_ikp, jsp, ik)
   !    use m_types
   !    use m_judft
   !    implicit none
   !    type(t_atoms), intent(in)    :: atoms
   !    type(t_cell), intent(in)     :: cell
   !    type(t_enpara), intent(in)   :: enpara
   !    type(t_input), intent(in)    :: input
   !    type(t_noco), intent(in)     :: noco
   !    type(t_hybinp), intent(in)   :: hybinp
   !    type(t_hybdat), intent(in)   :: hybdat
   !    type(t_kpts), intent(in)     :: kpts
   !    type(t_sym), intent(in)      :: sym
   !    type(t_oneD), intent(in)     :: oneD
   !    type(t_mat), intent(in)      :: zmat_ikp ! zmat of parent k-point
   !    integer, intent(in)          :: jsp
   !    integer, intent(in)          :: ik       ! k-point index
   !
   !
   !
   !    complex, allocatable :: acof(:, :, :), bcof(:, :, :), ccof(:, :, :, :)
   !    complex              :: cmt(input%neig, hybdat%maxlmindx, atoms%nat)
   !    complex              :: cmthlp(input%neig, hybdat%maxlmindx, atoms%nat)
   !
   !    integer :: ikp, nbands, ok(3), maxp ! index of parent k-point
   !
   !    type(t_lapw)  :: lapw_ik, lapw_ikp
   !    type(t_usdus) :: usdus
   !
   !    ikp = kpts%bkp(ik)
   !    nbands = hybdat%nbands(ikp)
   !
   !    allocate(acof(nbands, 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok(1))
   !    allocate(bcof(nbands, 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok(2))
   !    allocate(ccof(-atoms%llod:atoms%llod, nbands, atoms%nlod, atoms%nat), stat=ok(3))
   !    ! allocate(cmt(input%neig, hybdat%maxlmindx, atoms%nat), stat=ok)
   !    ! allocate(cmthlp(input%neig, hybdat%maxlmindx, atoms%nat), stat=ok)
   !
   !    CALL usdus%init(atoms, input%jspins)
   !    DO itype = 1, atoms%ntype
   !       DO l = 0, atoms%lmax(itype)
   !          CALL radfun(l, itype, jsp, enpara%el0(l,itype, jsp)), vr(:, itype, jsp), &
   !                    atoms, u(:, :, l), du(:, :, l), usdus, nodem, noded, wronk)
   !       END DO
   !
   !       IF (atoms%nlo(itype) >= 1) THEN
   !          CALL radflo(atoms, itype, jsp, ello_eig, vr(:, itype, jsp), u, du, mpi, usdus, uuilon, duilon, ulouilopn, flo)
   !       END IF
   !    END DO
   !
   !    if(any(ok /= 0)) then
   !       maxp = maxloc(abs(ok))
   !       call judft_error("Error in allocating array no" // int2str(maxp))
   !    endif
   !    acof = 0; bcof = 0; ccof = 0
   !
   !    CALL lapw_ik%init(input, noco, kpts, atoms, sym, ik, cell, sym%zrfs)
   !    CALL lapw_ikp%init(input, noco, kpts, atoms, sym, ikp, cell, sym%zrfs)
   !
   !    lapw_ikp%nmat = lapw_ikp%nv(jsp) + atoms%nlotot
   !
   !    CALL abcof(input, atoms, sym, cell, lapw_ikp, nbands, usdus, noco, jsp, &
   !               oneD, acof, bcof, ccof, zmat_ikp)
   !    CALL hyb_abcrot(hybinp, atoms, nbands, sym, acof, bcof, ccof)
   !
   ! end subroutine calc_cmt
end module m_calc_cmt
