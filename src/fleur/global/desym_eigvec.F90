!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_desym_eigvec
   use m_juDFT
   use m_constants
   use m_eig66_io
   use m_genMTBasis
   use m_trafo, only: waveftrafo_gen_zmat, waveftrafo_gen_cmt
   use m_types_abc
   use m_types_atoms
   use m_types_cell
   use m_types_desym_map
   use m_types_enpara
   use m_types_fleurinput
   use m_types_hybinp
   use m_types_input
   use m_types_lapw
   use m_types_mat
   use m_types_misc
   use m_types_mpdata
   use m_types_mpi
   use m_types_noco
   use m_types_nococonv
   use m_types_potden
   use m_types_radfun
   use m_types_sym
   use m_types_usdus
   implicit none
   private
   public :: desym_eigvec

contains

   subroutine desym_eigvec(fi, enpara, nococonv, vTot, eig_id, results, &
                           fmpi_nosym, fi_nosym, enpara_nosym, nococonv_nosym, vTot_nosym, map, &
                           eig_id_nosym, results_nosym)
      !! Rotates the IBZ eigenvectors of eig_id onto the full k-mesh of the setup without symmetry (eig_id_nosym).
      type(t_fleurinput), intent(in)  :: fi, fi_nosym
      type(t_enpara),     intent(in)  :: enpara, enpara_nosym
      type(t_nococonv),   intent(in)  :: nococonv, nococonv_nosym
      type(t_potden),     intent(in)  :: vTot, vTot_nosym
      integer,            intent(in)  :: eig_id
      type(t_results),    intent(in)  :: results
      type(t_mpi),        intent(in)  :: fmpi_nosym
      type(t_desym_map),  intent(in)  :: map
      integer,            intent(out) :: eig_id_nosym
      type(t_results),    intent(out) :: results_nosym

      type(t_hybinp) :: hybinp
      type(t_mpdata) :: mpdata
      type(t_usdus)  :: usdus, usdus_nosym
      type(t_lapw)   :: lapw_k, lapw_nosym
      type(t_mat)    :: zk, zk_c, z_nosym
      complex, allocatable :: cmt_k(:,:,:), cmt_rot(:,:,:), c_phase(:)
      integer, allocatable :: ev_list(:)
      integer :: jsp, ikk, ik, kp, iop, neig, nb, n, maxlmindx

      if (fi%noco%l_noco .or. fi%noco%l_soc) call juDFT_error("Eigenvector rotation is not implemented for noco/SOC.", calledby="desym_eigvec")

      call timestart("desym_eigvec")

      call hybinp%gen_map_all(fi%atoms, fi%sym, fi%input%film)
      call hybinp%gen_d_wgn2(fi%atoms, fi%cell, fi%sym)
      call mpdata%set_num_radfun_per_l(fi%atoms)
      maxlmindx = max_lmindx(fi%atoms)

      call desym_results(fi, fi_nosym, map, results, results_nosym)

      eig_id_nosym = open_eig(fmpi_nosym%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi_nosym%kpts%nkpt, fi%input%jspins, .false., &
                              .true., .false., .false., .false., .false., fmpi_nosym%n_size)

      do jsp = 1, fi%input%jspins
         call mt_basis(fi, enpara, vTot, fmpi_nosym, jsp, usdus)
         call mt_basis(fi_nosym, enpara_nosym, vTot_nosym, fmpi_nosym, jsp, usdus_nosym)

         do ikk = 1, size(fmpi_nosym%k_list)
            ik   = fmpi_nosym%k_list(ikk)
            kp   = map%kpt_parent(ik)
            iop  = map%kpt_op(ik)
            neig = results%neig(kp, jsp)
            ev_list = [(n, n = fmpi_nosym%n_rank + 1, neig, fmpi_nosym%n_size)]
            nb = size(ev_list)

            call lapw_k%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, kp, fi%cell)
            call lapw_nosym%init(fi_nosym%input, fi_nosym%noco, nococonv_nosym, fi_nosym%kpts, fi_nosym%atoms, fi_nosym%sym, ik, fi_nosym%cell)
            if (lapw_k%nv(jsp) /= lapw_nosym%nv(jsp)) call juDFT_error("Rotated k-point has a different number of plane waves.", calledby="desym_eigvec")
            lapw_k%nmat     = lapw_k%nv(jsp) + fi%atoms%nlotot
            lapw_nosym%nmat = lapw_nosym%nv(jsp) + fi_nosym%atoms%nlotot

            call z_nosym%init(.false., lapw_nosym%nmat, nb)

            if (nb > 0) then
               call zk%init(fi%input%l_real, lapw_k%nmat, nb)
               call read_eig(eig_id, kp, jsp, list=ev_list, zmat=zk)

               ! MT coefficients of the parent states in the global frame, rotated to k
               allocate(cmt_k(nb, maxlmindx, fi%atoms%nat), cmt_rot(nb, maxlmindx, fi%atoms%nat), c_phase(nb))
               call mt_coefficients(fi%input, fi%atoms, fi%sym, fi%cell, fi%noco, nococonv, usdus, lapw_k, zk, jsp, nb, cmt_k, hybinp)
               c_phase = cmplx_0
               call waveftrafo_gen_cmt(cmt_k, c_phase, .false., kp, iop, fi%atoms, mpdata, hybinp, fi%kpts, fi%sym, nb, cmt_rot)

               ! Plane-wave part; complex because the setup without symmetry has no inversion
               call zk_c%init(.false., zk%matsize1, nb)
               if (zk%l_real) then
                  zk_c%data_c = cmplx(zk%data_r, 0.0)
               else
                  zk_c%data_c = zk%data_c
               end if
               call waveftrafo_gen_zmat(zk_c, kp, iop, fi%kpts, fi%sym, jsp, nb, lapw_k, lapw_nosym, z_nosym)

               call lo_coefficients(fi_nosym, nococonv_nosym, usdus_nosym, lapw_nosym, jsp, nb, cmt_rot, z_nosym)

               deallocate(cmt_k, cmt_rot, c_phase)
               call zk%free()
               call zk_c%free()
            end if

            call write_eig(eig_id_nosym, ik, jsp, nb, neig, results%eig(:neig, kp, jsp), n_start=fmpi_nosym%n_size, &
                           n_end=fmpi_nosym%n_rank, zmat=z_nosym)
            call z_nosym%free()
         end do
      end do

      call timestop("desym_eigvec")

   end subroutine desym_eigvec

   subroutine desym_results(fi, fi_nosym, map, results, results_nosym)
      !! Copies results to the full k-mesh and the atom types of the setup without symmetry; w_iks is reweighted.
      type(t_fleurinput), intent(in)  :: fi, fi_nosym
      type(t_desym_map),  intent(in)  :: map
      type(t_results),    intent(in)  :: results
      type(t_results),    intent(out) :: results_nosym

      integer :: ik, kp
      real    :: wsym, wnosym

      results_nosym = results
      results_nosym%neig  = results%neig(map%kpt_parent, :)
      results_nosym%eig   = results%eig(:, map%kpt_parent, :)
      results_nosym%w_iks = results%w_iks(:, map%kpt_parent, :)
      do ik = 1, fi_nosym%kpts%nkpt
         kp = map%kpt_parent(ik)
         wsym   = fi%kpts%wtkpt(kp) / sum(fi%kpts%wtkpt)
         wnosym = fi_nosym%kpts%wtkpt(ik) / sum(fi_nosym%kpts%wtkpt)
         results_nosym%w_iks(:, ik, :) = results_nosym%w_iks(:, ik, :) * wnosym / wsym
      end do

      if (allocated(results%force))     results_nosym%force     = results%force(:, map%parent_type, :)
      if (allocated(results%force_old)) results_nosym%force_old = results%force_old(:, map%parent_type)
      if (allocated(results_nosym%unfolding_weights)) deallocate(results_nosym%unfolding_weights)
      if (allocated(results_nosym%w_iksRDMFT))        deallocate(results_nosym%w_iksRDMFT)
      ! The Wannier gauge lives on bkf, which is the k order without symmetry
      if (allocated(results_nosym%U_mat)) then
         if (size(results_nosym%U_mat, 3) /= fi_nosym%kpts%nkpt) deallocate(results_nosym%U_mat)
      end if
      if (allocated(results_nosym%U_dis)) then
         if (size(results_nosym%U_dis, 3) /= fi_nosym%kpts%nkpt) deallocate(results_nosym%U_dis)
      end if

   end subroutine desym_results

   subroutine mt_basis(fi, enpara, vTot, fmpi, jsp, usdus)
      !! Radial basis integrals (usdus) of all atom types for spin jsp.
      type(t_fleurinput), intent(in)  :: fi
      type(t_enpara),     intent(in)  :: enpara
      type(t_potden),     intent(in)  :: vTot
      type(t_mpi),        intent(in)  :: fmpi
      integer,            intent(in)  :: jsp
      type(t_usdus),      intent(out) :: usdus

      real, allocatable :: f(:,:,:), g(:,:,:), flo(:,:,:)
      integer :: itype

      allocate(f(fi%atoms%jmtd, 2, 0:fi%atoms%lmaxd), g(fi%atoms%jmtd, 2, 0:fi%atoms%lmaxd), flo(fi%atoms%jmtd, 2, fi%atoms%nlod))
      call usdus%init(fi%atoms, fi%input%jspins)
      do itype = 1, fi%atoms%ntype
         call genMTBasis(fi%atoms, enpara, vTot, fmpi, itype, jsp, usdus, f, g, flo, l_writeArg=.false.)
      end do

   end subroutine mt_basis

   subroutine mt_coefficients(input, atoms, sym, cell, noco, nococonv, usdus, lapw, zmat, jsp, nb, cmt, hybinp)
      !! MT coefficients cmt(band, lmindx, atom) of zmat in the hybrid-code layout (l, m, radial function).
      !! With hybinp, coefficients given in rotated local frames are transformed to the global frame.
      type(t_input),    intent(in)           :: input
      type(t_atoms),    intent(in)           :: atoms
      type(t_sym),      intent(in)           :: sym
      type(t_cell),     intent(in)           :: cell
      type(t_noco),     intent(in)           :: noco
      type(t_nococonv), intent(in)           :: nococonv
      type(t_usdus),    intent(in)           :: usdus
      type(t_lapw),     intent(in)           :: lapw
      type(t_mat),      intent(in)           :: zmat
      integer,          intent(in)           :: jsp, nb
      complex,          intent(out)          :: cmt(:,:,:)
      type(t_hybinp),   intent(in), optional :: hybinp

      type(t_abc)    :: abc
      type(t_radfun) :: rf
      integer :: itype, na, iatom, indx, l, m, i

      cmt = cmplx_0
      do itype = 1, atoms%ntype
         call abc%init(input, atoms, nb, itype)
         call rf%from_usdus(atoms, usdus, itype)
         call abc%calc_abc(input, atoms, sym, cell, lapw, nb, rf, noco, nococonv, jsp, itype, zmat)
         if (present(hybinp)) call abc%rot_to_unrotated(hybinp, atoms, sym, itype)
         do na = 1, atoms%neq(itype)
            iatom = atoms%firstAtom(itype) + na - 1
            indx = 0
            do l = 0, atoms%lmax(itype)
               do m = -l, l
                  do i = 1, abc%n_r(l)
                     indx = indx + 1
                     cmt(:, indx, iatom) = ImagUnit**l * abc%cof(:, l*(l+1)+m, i, na)
                  end do
               end do
            end do
         end do
      end do

   end subroutine mt_coefficients

   subroutine lo_coefficients(fi, nococonv, usdus, lapw, jsp, nb, cmt, zmat)
      !! Sets the local-orbital rows of zmat such that its LO coefficients in the MT spheres equal those of cmt.
      !! Only local orbitals contribute to the LO slot, so each atom and LO gives a (2l+1)x(2l+1) linear system.
      type(t_fleurinput), intent(in)    :: fi
      type(t_nococonv),   intent(in)    :: nococonv
      type(t_usdus),      intent(in)    :: usdus
      type(t_lapw),       intent(in)    :: lapw
      integer,            intent(in)    :: jsp, nb
      complex,            intent(in)    :: cmt(:,:,:)
      type(t_mat),        intent(inout) :: zmat

      type(t_abc)    :: abc
      type(t_radfun) :: rf
      type(t_mat)    :: zunit
      complex, allocatable :: amat(:,:), rhs(:,:)
      integer, allocatable :: ipiv(:)
      integer :: itype, iatom, lo, l, m, slot, ikv, c0, nv, n_r(0:fi%atoms%lmaxd), offset, info, i

      if (fi%atoms%nlotot == 0) return
      nv = lapw%nv(jsp)

      ! Unit vectors on the LO rows: calc_abc then gives the MT coefficients of each LO basis function
      call zunit%init(.false., zmat%matsize1, fi%atoms%nlotot)
      do i = 1, fi%atoms%nlotot
         zunit%data_c(nv + i, i) = cmplx_1
      end do

      do itype = 1, fi%atoms%ntype
         if (fi%atoms%nlo(itype) == 0) cycle
         iatom = fi%atoms%firstAtom(itype)
         call abc%init(fi%input, fi%atoms, fi%atoms%nlotot, itype)
         call rf%from_usdus(fi%atoms, usdus, itype)
         call abc%calc_abc(fi%input, fi%atoms, fi%sym, fi%cell, lapw, fi%atoms%nlotot, rf, fi%noco, nococonv, jsp, itype, zunit)
         n_r = fi%atoms%num_radial_functions_per_l(itype)

         do lo = 1, fi%atoms%nlo(itype)
            l    = fi%atoms%llo(lo, itype)
            slot = fi%atoms%slot_of_lo(lo, itype)
            c0   = lapw%index_lo(lo, iatom)
            if (lapw%nkvec(lo, iatom) /= 2*l + 1) call juDFT_error("Unexpected number of LO basis functions.", calledby="lo_coefficients")
            offset = sum([((2*i + 1)*n_r(i), i = 0, l - 1)])

            allocate(amat(2*l + 1, 2*l + 1), rhs(2*l + 1, nb), ipiv(2*l + 1))
            do m = -l, l
               do ikv = 1, 2*l + 1
                  amat(m + l + 1, ikv) = abc%cof(c0 + ikv, l*(l+1)+m, slot, 1)
               end do
               rhs(m + l + 1, :) = cmt(:, offset + (m + l)*n_r(l) + slot, iatom) / ImagUnit**l
            end do
            call zgesv(2*l + 1, nb, amat, 2*l + 1, ipiv, rhs, 2*l + 1, info)
            if (info /= 0) call juDFT_error("LO coefficient system is singular.", calledby="lo_coefficients")
            zmat%data_c(nv + c0 + 1:nv + c0 + 2*l + 1, :) = rhs
            deallocate(amat, rhs, ipiv)
         end do
      end do

      call zunit%free()

   end subroutine lo_coefficients

   integer function max_lmindx(atoms)
      !! Largest number of (l, m, radial function) entries of any atom type.
      type(t_atoms), intent(in) :: atoms

      integer :: itype, l, n_r(0:atoms%lmaxd)

      max_lmindx = 0
      do itype = 1, atoms%ntype
         n_r = atoms%num_radial_functions_per_l(itype)
         max_lmindx = max(max_lmindx, sum([((2*l + 1)*n_r(l), l = 0, atoms%lmax(itype))]))
      end do

   end function max_lmindx

end module m_desym_eigvec
