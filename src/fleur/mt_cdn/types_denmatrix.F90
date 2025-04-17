!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_types_denmatrix
   use m_judft
   implicit none
   type:: t_denmatrix
      integer, private:: itype = 0
      logical         :: l_triang = .false.
      complex, allocatable  :: mat(:, :, :, :, :) !!radial function (index1,index2,l1,l2,lh)
   contains
      procedure, pass :: init
      procedure, pass :: rhonmt
      procedure, pass :: l_like_charge
      procedure, pass :: to_full_density
      procedure, pass  :: mpi_collect
   end type

contains

   subroutine mpi_collect(this)
      implicit none
      class(t_denmatrix), intent(inout):: this
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%mat, size(this%mat), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine
   pure subroutine init(this, itype, atoms, input, sphhar)
      use m_types
      use m_types_radfun
      implicit none
      class(t_denmatrix), intent(inout):: this
      integer, intent(in):: itype
      type(t_atoms), intent(in):: atoms
      type(t_input), intent(in):: input
      type(t_sphhar), intent(in):: sphhar

      type(t_radfun):: radfun
      integer:: max_rad_funct

      this%itype = itype

      call radfun%init(atoms, input, itype)
      max_rad_funct = maxval(radfun%n_r)
      if (allocated(this%mat)) deallocate (this%mat)
      Allocate (this%mat(max_rad_funct, max_rad_funct, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:sphhar%nlhd))
      this%mat = 0.0
   end subroutine

   function l_like_charge(this, radfun, ispin, lmax) result(qmtl)
      use m_types_radfun
      use m_constants
      implicit none
      class(t_denmatrix), intent(in)   :: this
      type(t_radfun), intent(in)       :: radfun
      integer, intent(in)              :: ispin, lmax
      real :: qmtl(0:lmax)

      integer:: l, i, j
      qmtl = 0.0
      !l-decomposed density
      do l = 0, lmax
         DO i = 1, radfun%n_r(l)
            DO j = 1, radfun%n_r(l)
               qmtl(l) = qmtl(l) + real(this%mat(i, j, l, l, 0))*radfun%integral(i, j, l, ispin, ispin)*sfp_const
            end do
         END DO
      END DO
   end function

   subroutine rhonmt(this, atoms, sphhar, we, ne, itype, sym, abc, abc1, abc1m)
    !! Subroutine to construct all non-spherical MT density coefficients (for a
    !! density perturbation) without LOs in one routine. The spin input dictates,
    !! which element is gonna be built.
    !! The coefficients are of the form:
    !! $$\begin{aligned}
    !! d_{l',l,L,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha} &= \sum_{\nu\boldsymbol{k}}\sum_{m\mu(L)}\\
    !! &*  c_{L,\mu}^{*}G_{l,l''(L),l'}^{m,m''(\mu),m-m''(\mu)}A_{l',m-m''(\mu),\lambda'}^{\sigma_{\alpha}',\nu\boldsymbol{k}*}\\
    !! &* (2\tilde{f}_{\nu\boldsymbol{k}}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}\boldsymbol{q},j,\beta~(1)}+\tilde{f}_{\nu\boldsymbol{k}\boldsymbol{q}}^{j,\beta~(1)}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}})
    !! \end{aligned}$$
    !! The k-point loop is performed outside this routine. In contrast to older
    !! routines, the arrays uunmt etc. and uunmt21 etc. are merged into one
    !! spin and u-order dependent array.
    !!
    !! \(G_{l,l'',l'}^{m,m'',m'}\): Gaunt coefficient
    !!
    !! \(\sigma_{\alpha}(')\): local spin indices \(\rightarrow\) ilSpinPr, ilSpin
    !!
    !! \(\lambda(')\): order of the radial basis function (0: u, 1: d)
    !!
    !! \(L\): Lattice harmonic index
    !!
    !! \(\mu(L)\): Lattice harmonic member index
    !!
    !! \(c_{L,\mu}\): Lattice harmonic coefficients
    !!
    !! \(\nu\boldsymbol{k}\): State index (k-point and number of state)
    !!
    !! \(\boldsymbol{q},j,\beta\): Perturbation index; q-point, direction and atom
    !!
    !! \(\tilde{f}_{\nu\boldsymbol{k}}\): (Smeared) occupation number [perturbed for \(\tilde{f}^{(1)}\)]
    !!
    !! \(A\): Summed matching coefficients and eigenvectors [perturbed for \(A^{(1)}\)]
      use m_types
      use m_types_abc
      use m_gaunt
      implicit none
      class(t_denmatrix)         :: this
      type(t_sym), intent(IN)    :: sym
      type(t_sphhar), intent(IN)    :: sphhar
      type(t_atoms), intent(IN)    :: atoms

      integer, intent(IN)    :: ne, itype
      real, intent(IN)    :: we(ne)    !! \(\tilde{f}_{\nu\boldsymbol{k}}\)
      !real, intent(IN)    :: we1(ne)   !! \(\tilde{f}_{\nu\boldsymbol{k}\boldsymbol{q}}^{j,\beta~(1)}\)
      !logical, intent(IN)    :: l_dfpt

      type(t_abc), intent(IN)    :: abc  !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}}\)
      type(t_abc), intent(IN)    :: abc1 !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}\boldsymbol{q},j,\beta~(1)}\)

      type(t_abc), optional, intent(IN) :: abc1m !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}\boldsymbol{-q},j,\beta~(1)}\)

      complex :: coef, cil, cmv
      complex :: temp(ne)
      real:: fact
      integer :: jmem, l, lh, llp, llpmax, lm, lmp, lp, lv, m, mp, mv, na, natom, nn, ns, nt, lphi, lplow, icoef, jcoef, mpp,lpp,lpmin0,lpmax0,lpmin,lo,lop,lpmax
      integer :: n_l(atoms%nlod)
      !TODO these are needed for DFPT!?
      logical :: l_minusq, l_less_effort, l_gamma, l_dfpt
      real, allocatable:: we1(:)!(nobd)
      l_minusq = present(abc1m)
      l_less_effort = .true. !TODO In case of LOs this does not work
      l_gamma = .false.
      l_dfpt = .false.
      this%l_triang = l_less_effort
      ns = sym%ntypsy(atoms%firstAtom(itype))
      do lh = 0, sphhar%nlh(ns)
         lv = sphhar%llh(lh, ns)
         do jmem = 1, sphhar%nmem(lh, ns)
            mv = sphhar%mlh(jmem, lh, ns)
            cmv = conjg(sphhar%clnu(jmem, lh, ns))
            do l = 0, atoms%lmax(itype)
               m_loop: do m = -l, l
                  lm = l*(l + 1) + m
                  mp = m - mv
                  !maximum value of lp
                  lphi = l + lv
                  !---> check that lphi is smaller than the max l of the
                  !---> wavefunction expansion
                  lphi = min(lphi, atoms%lmax(itype))
                  !--->  make sure that l + l'' + lphi is even
                  lphi = lphi - mod(l + lv + lphi, 2)
                  if (this%l_triang) lphi = l - mod(lv, 2)
                  lplow = abs(l - lv)
                  lplow = max(lplow, abs(mp))
                  !---> make sure that l + l'' + lplow is even
                  lplow = lplow + mod(abs(lphi - lplow), 2)
                  if (lplow .gt. lphi) cycle m_loop
                  do lp = lplow, lphi, 2
                     cil = ImagUnit**(l - lp)
                     lmp = lp*(lp + 1) + mp
                     if (lmp > lm .and. this%l_triang) cycle m_loop
                     coef = cmv*cil*gaunt1(l, lv, lp, m, mv, mp, atoms%lmaxd)
                     !IF (ABS(coef) .LT. 1e-12 ) CYCLE
                     do nt = 1, atoms%neq(itype)
                        ! uu/du
                        do icoef = 1, size(abc%cof, 3) !Loop over radial functions/cofs (usually 0=u and 1=\dot u)
                           temp(:) = coef*we(:)*abc1%cof(:, lm, icoef, nt) ! If not DFPT, this is the base case for rhonmt(21)
                           if (lmp /= lm .and. this%l_triang) temp(:) = temp(:)*2.0
                           if (l_dfpt) then
                              if (.not. l_minusq) temp(:) = temp(:)*2.0
                              if (l_gamma) temp(:) = temp(:) + coef*we1(:)*abc%cof(:, lm, 0, nt)
                           end if
                           do jcoef = 1, size(abc%cof, 3) !Loop over radial functions/cofs (usually 0=u and 1=\dot u)
                              this%mat(jcoef, icoef, l, lp, lh) = this%mat(jcoef, icoef, l, lp, lh) &
                                                    & + dot_product(abc%cof(:ne, lmp, jcoef, nt), temp(:ne))
                              if (l_minusq) then
                                 this%mat(jcoef, icoef, l, lp, lh) = this%mat(jcoef, icoef, l, lp, lh) &
                                                          & + dot_product(abc1m%cof(:ne, lmp, jcoef, nt), &
                                                               & abc%cof(:, lm, icoef, nt))
                              end if
                           end do
                        end do

                     end do ! nt

                  end do!lp
               end do m_loop ! m
            end do ! jmem
         end do ! l
      end do ! lh
   end subroutine rhonmt

   subroutine to_full_density(denmat, ispin, ispinpr, itype, input, sphhar, atoms, noco, sym, radfun, rho, rhoIm, moments)
      !! Current situation:
      !!
      !! This routine calculates density contributions
      !! $$\rho_{L}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha}(r)=
      !! \sum_{l',l,\lambda',\lambda,s}d_{l',l,L,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha}
      !! u_{l',\lambda',s}^{\sigma_{\alpha}',\alpha}(r)u_{l,\lambda,s}^{\sigma_{\alpha},\alpha}(r)$$
      !! \(s\) is the index for the big/small components yielded by the
      !! scalar-relativistic Schrödinger equation.

      use m_types
      use m_types_radfun

      implicit none
      CLASS(t_denmatrix), intent(IN)  :: denmat
      integer, intent(in)             :: itype, ispin, ispinpr
      type(t_input), intent(IN)      :: input
      type(t_sphhar), intent(IN)     :: sphhar
      type(t_atoms), intent(IN)      :: atoms
      type(t_sym), intent(IN)        :: sym
      type(t_noco), intent(IN)       :: noco
      type(t_radfun), intent(in)      :: radfun
      real, intent(inout)             :: rho(:, 0:, :, :)

      real, intent(inout), optional    :: rhoIm(:, 0:)
      type(t_moments), optional, intent(INOUT) :: moments

      integer :: lh, l, lp, llp, j, i, ii, spin
      complex :: cs
      integer, parameter:: lcf = 3

      spin = merge(ispin, 3, ispin == ispinpr)
      if (spin == 3 .and. .not. noco%l_mperp) return !no local off-diagonal pert
      do lh = 0, sphhar%nlh(sym%ntypsy(atoms%firstAtom(itype)))
         do l = 0, atoms%lmax(itype)
            do lp = 0, merge(l, atoms%lmax(itype), denmat%l_triang)
               llp = (l*(l + 1))/2 + lp

               !if (lh > 0 .and. atoms%l_outputCFpot(itype) .and. atoms%l_outputCFremove4f(itype) &
               !    .and. (l == lcf .and. lp == lcf)) cycle !Exclude non-spherical contributions for CF

               do j = 1, atoms%jri(itype)
                  cs = 0.0
                  do i = 1, radfun%n_r(l) !Loop over radial functions
                     do ii = 1, radfun%n_r(lp)
                        cs = cs + denmat%mat(ii, i, l, lp, lh) &   !density matrix
                             *(radfun%r(i, j, 1, l, ispinpr)*radfun%R(ii, j, 1, lp, ispin) &  !large components
                               + radfun%R(i, j, 2, l, ispinpr)*radfun%R(ii, j, 2, lp, ispin))    !small components
                     end do
                  end do
                  rho(j, lh, itype, spin) = rho(j, lh, itype, spin) + real(cs)/atoms%neq(itype)
                  if (spin == 3) rho(j, lh, itype, 4) = rho(j, lh, itype, 4) + aimag(cs)/atoms%neq(itype) !Store imaginary part as 4th spin
                  if (present(rhoIm)) rhoIm(j, lh) = rhoIm(j, lh) + aimag(cs)/atoms%neq(itype)
                  if ((l <= input%lResMax) .and. (lp <= input%lResMax) .and. ispin == ispinpr .and. present(moments)) &
                   moments%rhoLRes(j, lh, llp, itype, ispin) = moments%rhoLRes(j, lh, llp, itype, ispin) + real(cs)/atoms%neq(itype)

               end do
            end do
         end do

      END DO
   end subroutine to_full_density
end module m_types_denmatrix
