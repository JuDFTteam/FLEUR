module m_calc_cmt

contains
   subroutine calc_cmt(atoms, cell, input, noco, nococonv, hybinp, hybdat, mpdata, kpts, &
                       sym,   zmat_ikp, jsp, ik, c_phase, cmt_out, submpi)
      use m_types
      use m_judft
      USE m_hyb_abcrot
      USE m_abcof
      use m_constants
      use m_trafo, only: waveftrafo_gen_cmt
      use m_io_hybrid
      use m_divide_most_evenly 
      implicit none
      type(t_atoms), intent(in)       :: atoms
      type(t_cell), intent(in)        :: cell
      type(t_input), intent(in)       :: input
      type(t_noco), intent(in)        :: noco
      type(t_nococonv), intent(in)    :: nococonv
      type(t_hybinp), intent(in)      :: hybinp
      type(t_hybdat), intent(in)      :: hybdat
      type(t_mpdata), intent(in)      :: mpdata
      type(t_kpts), intent(in)        :: kpts
      type(t_sym), intent(in)         :: sym
       
      type(t_mat), intent(in), target :: zmat_ikp ! zmat of parent k-point
                                               ! not sure wether zmat works aswell
      integer, intent(in)          :: jsp
      integer, intent(in)          :: ik       ! k-point
      complex, intent(in)          :: c_phase(:)

      complex, intent(inout)       :: cmt_out(:,:,:)
      type(t_hybmpi), intent(in), optional :: submpi
      complex, allocatable :: acof(:,:,:), bcof(:,:,:), ccof(:,:,:,:)
      complex, allocatable :: cmt(:,:,:)

      integer :: ikp, nbands, ok(4) ! index of parent k-point
      integer :: iatom, itype, indx, i, j, idum, iop, l, ll, lm, m
      integer :: map_lo(atoms%nlod)
      integer, allocatable :: start_idx(:), psize(:)
      integer :: my_psz, my_start, ierr

      complex :: cdum
      type(t_lapw)  :: lapw_ik, lapw_ikp
      type(t_mat), target   :: tmp
      type(t_mat), pointer  :: mat_ptr

      call timestart("calc_cmt")
      ikp = kpts%bkp(ik)
      nbands = zmat_ikp%matsize2

      if(present(submpi)) then
         call divide_most_evenly(nbands, submpi%size, start_idx, psize)
         my_psz = psize(submpi%rank+1)
         my_start = start_idx(submpi%rank+1)
      else
         my_psz = nbands
         my_start = 1
      endif

      call timestart("alloc abccof")
      allocate(acof(my_psz, 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok(1), source=cmplx_0)
      allocate(bcof(my_psz, 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok(2), source=cmplx_0)
      allocate(ccof(-atoms%llod:atoms%llod, my_psz, atoms%nlod, atoms%nat), stat=ok(3), source=cmplx_0)
      allocate(cmt(nbands, hybdat%maxlmindx, atoms%nat), stat=ok(4), source=cmplx_0)
      if(any(ok /= 0)) then
         call judft_error("Error in allocating abcof arrays")
      endif
      call timestop("alloc abccof")

      CALL lapw_ik%init(input, noco, nococonv, kpts, atoms, sym, ik, cell)
      CALL lapw_ikp%init(input, noco, nococonv, kpts, atoms, sym, ikp, cell)

      lapw_ikp%nmat = lapw_ikp%nv(jsp) + atoms%nlotot
      if(my_psz /= nbands) then
         call tmp%init(zmat_ikp%l_real, zmat_ikp%matsize1, my_psz)
         if(tmp%l_real) then 
            tmp%data_r = zmat_ikp%data_r(:,my_start:my_start+my_psz-1)
         else
            tmp%data_c = zmat_ikp%data_c(:,my_start:my_start+my_psz-1)
         endif
         mat_ptr => tmp 
      else 
         mat_ptr => zmat_ikp
      endif

      CALL abcof(input, atoms, sym, cell, lapw_ikp, my_psz, hybdat%usdus, noco, nococonv,&
                 jsp,   acof, bcof, ccof, mat_ptr)
      CALL hyb_abcrot(hybinp, atoms, my_psz, sym, acof, bcof, ccof)

      call timestart("copy to cmt")
      !$OMP parallel do default(none) private(iatom, itype, indx, l, ll, cdum, idum, map_lo, j, m, lm, i) &
      !$OMP shared(atoms, mpdata, cmt, acof, bcof, ccof, my_start, my_psz)
      do iatom = 1,atoms%nat 
         itype = atoms%itype(iatom)
         indx = 0
         DO l = 0, atoms%lmax(itype)
            ll = l*(l + 1)
            cdum = ImagUnit**l

            ! determine number of local orbitals with quantum number l
            ! map returns the number of the local orbital of quantum
            ! number l in the list of all local orbitals of the atom type
            idum = 0
            map_lo = 0
            IF (mpdata%num_radfun_per_l(l, itype) > 2) THEN
               DO j = 1, atoms%nlo(itype)
                  IF (atoms%llo(j, itype) == l) THEN
                     idum = idum + 1
                     map_lo(idum) = j
                  END IF
               END DO
            END IF

            DO M = -l, l
               lm = ll + M
               DO i = 1, mpdata%num_radfun_per_l(l, itype)
                  indx = indx + 1
                  IF (i == 1) THEN
                     cmt(my_start:my_start+my_psz-1, indx, iatom) = cdum*acof(:, lm, iatom)
                  ELSE IF (i == 2) THEN
                     cmt(my_start:my_start+my_psz-1, indx, iatom) = cdum*bcof(:, lm, iatom)
                  ELSE
                     idum = i - 2
                     cmt(my_start:my_start+my_psz-1, indx, iatom) = cdum*ccof(M, :, map_lo(idum), iatom)
                  END IF
               END DO
            END DO
         END DO
      END DO
      !$OMP end parallel do
      call timestop("copy to cmt")

#ifdef CPP_MPI
      if(my_psz /= nbands) then
         call timestart("allreduce cmt")
         do i = 1,atoms%nat
            call MPI_Allreduce(MPI_IN_PLACE, cmt(1,1,i), nbands*hybdat%maxlmindx, MPI_DOUBLE_COMPLEX, MPI_SUM, submpi%comm, ierr)
         enddo
         call timestop("allreduce cmt")
      endif
#endif

      ! write cmt at irreducible k-points in direct-access file cmt
      if(ik <= kpts%nkpt) then
         call timestart("copy out cmt")
         call zcopy(size(cmt), cmt, 1, cmt_out, 1)
         call timestop("copy out cmt")
      else
         iop = kpts%bksym(ik)

         call waveftrafo_gen_cmt(cmt, c_phase, zmat_ikp%l_real, ikp, iop, atoms, &
                                  mpdata, hybinp, kpts, sym, nbands, cmt_out)
      endif
      call timestop("calc_cmt")
   end subroutine calc_cmt
end module m_calc_cmt
