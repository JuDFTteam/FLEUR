module m_calc_cmt

contains
   subroutine calc_cmt(atoms, cell, input, noco, nococonv, hybinp, hybdat, mpdata, kpts, &
                       sym, oneD, zmat_ikp, jsp, ik, c_phase, cmt_out)
      use m_types
      use m_judft
      USE m_hyb_abcrot
      USE m_abcof
      use m_constants, only: ImagUnit
      use m_trafo, only: waveftrafo_gen_cmt
      use m_io_hybinp
      implicit none
      type(t_atoms), intent(in)    :: atoms
      type(t_cell), intent(in)     :: cell
      type(t_input), intent(in)    :: input
      type(t_noco), intent(in)     :: noco
      type(t_nococonv), intent(in) :: nococonv
      type(t_hybinp), intent(in)   :: hybinp
      type(t_hybdat), intent(in)   :: hybdat
      type(t_mpdata), intent(in)   :: mpdata
      type(t_kpts), intent(in)     :: kpts
      type(t_sym), intent(in)      :: sym
      type(t_oneD), intent(in)     :: oneD
      type(t_mat), intent(in)      :: zmat_ikp ! zmat of parent k-point
                                               ! not sure wether zmat works aswell
      integer, intent(in)          :: jsp
      integer, intent(in)          :: ik       ! k-point
      complex, intent(in)          :: c_phase(:)

      complex, intent(inout)       :: cmt_out(:,:,:)

      complex, allocatable :: acof(:,:,:), bcof(:,:,:), ccof(:,:,:,:)
      complex, allocatable :: cmt(:,:,:)

      integer :: ikp, nbands, ok(4) ! index of parent k-point
      integer :: iatom, itype, ieq, indx, i, j, idum, iop, l, ll, lm, m
      integer :: map_lo(atoms%nlod)

      complex :: cdum
      type(t_lapw)  :: lapw_ik, lapw_ikp

      call timestart("calc_cmt")

      ikp = kpts%bkp(ik)
      nbands = zmat_ikp%matsize2

      allocate(acof(nbands, 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok(1))
      allocate(bcof(nbands, 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), stat=ok(2))
      allocate(ccof(-atoms%llod:atoms%llod, nbands, atoms%nlod, atoms%nat), stat=ok(3))
      allocate(cmt(nbands, hybdat%maxlmindx, atoms%nat), stat=ok(4))

      if(any(ok /= 0)) then
         call judft_error("Error in allocating abcof arrays")
      endif
      acof = 0; bcof = 0; ccof = 0

      CALL lapw_ik%init(input, noco, nococonv, kpts, atoms, sym, ik, cell, sym%zrfs)
      CALL lapw_ikp%init(input, noco, nococonv, kpts, atoms, sym, ikp, cell, sym%zrfs)

      lapw_ikp%nmat = lapw_ikp%nv(jsp) + atoms%nlotot

      CALL abcof(input, atoms, sym, cell, lapw_ikp, nbands, hybdat%usdus, noco, nococonv,&
                 jsp, oneD, acof, bcof, ccof, zmat_ikp)
      CALL hyb_abcrot(hybinp, atoms, nbands, sym, acof, bcof, ccof)

      cmt = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
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
                        cmt(:, indx, iatom) = cdum*acof(:, lm, iatom)
                     ELSE IF (i == 2) THEN
                        cmt(:, indx, iatom) = cdum*bcof(:, lm, iatom)
                     ELSE
                        idum = i - 2
                        cmt(:, indx, iatom) = cdum*ccof(M, :, map_lo(idum), iatom)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      ! write cmt at irreducible k-points in direct-access file cmt
      if(ik <= kpts%nkpt) then
         cmt_out = cmt
      else
         iop = kpts%bksym(ik)

         call waveftrafo_gen_cmt(cmt, c_phase, zmat_ikp%l_real, ikp, iop, atoms, &
                                  mpdata, hybinp, kpts, sym, nbands, cmt_out)
      endif
      call timestop("calc_cmt")
   end subroutine calc_cmt
end module m_calc_cmt
