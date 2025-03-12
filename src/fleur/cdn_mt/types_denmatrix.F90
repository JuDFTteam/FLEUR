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
      complex,allocatable  :: mat(:, :, :, :, :) !!radial function (index1,index2,l1,l2,lh)
   contains
      procedure, pass :: init
      procedure, pass :: rhonmt
      procedure,pass  :: mpi_collect
   end type

contains

subroutine mpi_collect(this)
   implicit none
   class(t_denmatrix),intent(inout):: this
#ifdef CPP_MPI   
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,this%mat,size(this%mat),MPI_DOUBLE_COMPLEX, MPI_SUM,MPI_COMM_WORLD,ierr)
#endif   
end subroutine       
subroutine init(this,itype,atoms,input,sphhar)
   use m_types
   use m_types_radfun
   implicit none
   class(t_denmatrix):: this
   integer,intent(in):: itype
   type(t_atoms),intent(in):: atoms
   type(t_input),intent(in):: input
   type(t_sphhar),intent(in):: sphhar
   

   type(t_radfun):: radfun
   integer:: max_rad_funct

   this%itype=itype

   call radfun%init(atoms, input, itype)
   max_rad_funct=maxval(radfun%n_r)

   Allocate(this%mat(max_rad_funct,max_rad_funct,atoms%lmaxd,atoms%lmaxd,sphhar%nlhd))
   this%mat=0.0
end subroutine


subroutine rhonmt(this, atoms, sphhar, we, ne, itype, ilSpinPr, ilSpin, sym, eigVecCoeffs, eigVecCoeffs1,  eigVecCoeffs1m)
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
      use m_gaunt
      implicit none
      class(t_denmatrix)         :: this
      type(t_sym), intent(IN)    :: sym
      type(t_sphhar), intent(IN)    :: sphhar
      type(t_atoms), intent(IN)    :: atoms

      integer, intent(IN)    :: ne, itype
      integer, intent(IN)    :: ilSpinPr  !! \(\sigma_{\alpha}^{'}\)
      integer, intent(IN)    :: ilSpin    !! \(\sigma_{\alpha}\)
      real, intent(IN)    :: we(ne)    !! \(\tilde{f}_{\nu\boldsymbol{k}}\)
      !real, intent(IN)    :: we1(ne)   !! \(\tilde{f}_{\nu\boldsymbol{k}\boldsymbol{q}}^{j,\beta~(1)}\)
      !logical, intent(IN)    :: l_dfpt

      type(t_eigVecCoeffs), intent(IN)    :: eigVecCoeffs  !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}}\)
      type(t_eigVecCoeffs), intent(IN)    :: eigVecCoeffs1 !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}\boldsymbol{q},j,\beta~(1)}\)

    
      type(t_eigVecCoeffs), optional, intent(IN) :: eigVecCoeffs1m !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}\boldsymbol{-q},j,\beta~(1)}\)

      complex :: coef, cil, cmv
      complex :: temp(ne)
      real:: fact
      integer :: jmem, l, lh, llp, llpmax, lm, lmp, lp, lv, m, mp, mv, na, natom, nn, ns, nt, lphi, lplow, icoef, jcoef, mpp,lpp,lpmin0,lpmax0,lpmin,lo,lop,lpmax

      !TODO these are needed for DFPT!?
      logical :: l_minusq, l_less_effort, l_gamma, l_dfpt
      real,allocatable:: we1(:)!(nobd)
      l_minusq = present(eigVecCoeffs1m)
      l_less_effort = .true. !todo
      l_gamma = .false.
      l_dfpt = .false.

      ns = sym%ntypsy(atoms%firstAtom(itype))
      do lh = 0, sphhar%nlh(ns)
         lv = sphhar%llh(lh, ns)
         do jmem = 1, sphhar%nmem(lh, ns)
            mv = sphhar%mlh(jmem, lh, ns)
            cmv = conjg(sphhar%clnu(jmem, lh, ns))
            do l = 0, atoms%lmaxd
               m_loop: do m = -l, l
                  lm = l*(l + 1) + m
                  mp = m - mv
                  !maximum value of lp
                  lphi = l + lv
                  !---> check that lphi is smaller than the max l of the
                  !---> wavefunction expansion
                  lphi = min(lphi, atoms%lmaxd)
                  !--->  make sure that l + l'' + lphi is even
                  lphi = lphi - mod(l + lv + lphi, 2)
                  if (l_less_effort) lphi = l - mod(lv, 2)
                  lplow = abs(l - lv)
                  lplow = max(lplow, abs(mp))
                  !---> make sure that l + l'' + lplow is even
                  lplow = lplow + mod(abs(lphi - lplow), 2)
                  if (lplow .gt. lphi) cycle m_loop
                  do lp = lplow, lphi, 2
                     cil = ImagUnit**(l - lp)
                     lmp = lp*(lp + 1) + mp
                     if (lmp > lm .and. l_less_effort) cycle m_loop
                     coef = cmv*cil*gaunt1(l, lv, lp, m, mv, mp, atoms%lmaxd)
                     !IF (ABS(coef) .LT. 1e-12 ) CYCLE
                     do nt = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype)
                        ! uu/du
                        do icoef = lbound(eigVecCoeffs1%abcof, 3), ubound(eigVecCoeffs1%abcof, 3) !Loop over radial functions/abcofs (usually 0=u and 1=\dot u)
                           temp(:) = coef*we(:)*eigVecCoeffs1%abcof(:, lm, icoef, nt, ilSpin) ! If not DFPT, this is the base case for rhonmt(21)
                           if (lmp /= lm .and. l_less_effort) temp(:) = temp(:)*2.0
                           if (l_dfpt) then
                              if (.not. l_minusq) temp(:) = temp(:)*2.0
                              if (l_gamma) temp(:) = temp(:) + coef*we1(:)*eigVecCoeffs%abcof(:, lm, 0, nt, ilSpin)
                           end if
                           do jcoef = lbound(eigVecCoeffs1%abcof, 3), ubound(eigVecCoeffs1%abcof, 3) !Loop over radial functions/abcofs (usually 0=u and 1=\dot u)
                              this%mat(jcoef, icoef, l, lp, lh) = this%mat(jcoef, icoef, l, lp, lh) &
                                                    & + dot_product(eigVecCoeffs%abcof(:ne, lmp, jcoef, nt, ilSpinPr), temp(:ne))
                              if (l_minusq) then
                                 this%mat(jcoef, icoef, l, lp, lh) = this%mat(jcoef, icoef, l, lp, lh) &
                                                          & + dot_product(eigVecCoeffs1m%abcof(:ne, lmp, jcoef, nt, ilSpinPr), &
                                                               & eigVecCoeffs%abcof(:, lm, icoef, nt, ilSpin))
                              end if
                           end do
                        end do

                     end do ! nt

                  end do!lp
               end do m_loop ! m
            end do ! jmem
         end do ! l
      end do ! lh          
      !Now the LO-part
      ns = sym%ntypsy(atoms%firstatom(itype))
      do lh = 1, sphhar%nlh(ns)
         lpp = sphhar%llh(lh, ns)
         do jmem = 1, sphhar%nmem(lh, ns)
            mpp = sphhar%mlh(jmem, lh, ns)
            cmv = conjg(sphhar%clnu(jmem, lh, ns))
            do lo = 1, atoms%nlo(itype)
               l = atoms%llo(lo, itype)
               lpmin0 = abs(l - lpp)
               lpmax0 = l + lpp

               lpmax = min(lpmax0, atoms%lmax(itype))
               lpmax = lpmax - mod(l + lpp + lpmax, 2)
               do m = -l, l
                  mp = m - mpp
                  lpmin = max(lpmin0, abs(mp))
                  lpmin = lpmin + mod(l + lpp + lpmin, 2)
                  do lp = lpmin, lpmax, 2
                     lmp = lp*(lp + 1) + mp
                     fact = cmv*(ImagUnit**(l - lp))*gaunt1(l, lp, lpp, m, mp, mpp, atoms%lmaxd)
                     do na = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype)

                        temp = fact*we(:)*eigVecCoeffs1%ccof(m, :, lo, na, ilSpin)! If not DFPT, this is the base case for rhonmtlo
                        if (l_dfpt) then
                           temp = temp*2.0
                           if (l_gamma) then
                              temp = temp + fact*we1(:)*eigVecCoeffs%ccof(m, :, lo, na, ilSpin)
                           end if
                        end if
                        this%mat(2 + lo, 1, l, lp, lh) = this%mat(2 + lo, 1, l, lp, lh) &
                                                         + dot_product(eigVecCoeffs%abcof(:, lmp, 0, na, ilSpinPr), temp)
                        this%mat(2 + lo, 2, l, lp, lh) = this%mat(2 + lo, 2, l, lp, lh) &
                                               & + dot_product(eigVecCoeffs%abcof(:, lmp, 1, na, ilSpinPr), temp)

                     end do
                  end do

                  mp = m + mpp
                  lpmin = max(lpmin0, abs(mp))
                  lpmin = lpmin + mod(l + lpp + lpmin, 2)
                  do lp = lpmin, lpmax, 2
                     lmp = lp*(lp + 1) + mp
                     fact = cmv*(ImagUnit**(lp - l))*gaunt1(lp, l, lpp, mp, m, mpp, atoms%lmaxd)
                     do na = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype)

                        temp = fact*we(:)*eigVecCoeffs1%abcof(:, lmp, 0, na, ilSpin)! If not DFPT, this is the base case for rhonmtlo
                        if (l_dfpt) then
                           temp = temp*2.0
                           if (l_gamma) temp = temp + fact*we1(:)*eigVecCoeffs%abcof(:, lmp, 0, na, ilSpin)
                        end if
                        this%mat(1, 2 + lo, lp, l, lh) = this%mat(1, 2 + lo, lp, l, lh) &
                                          & + dot_product(eigVecCoeffs%ccof(m, :, lo, na, ilSpinPr), temp)
                        temp = fact*we(:)*eigVecCoeffs1%abcof(:, lmp, 1, na, ilSpin)! If not DFPT, this is the base case for rhonmtlo
                        if (l_dfpt) then
                           temp = temp*2.0
                           if (l_gamma) temp = temp + fact*we1(:)*eigVecCoeffs%abcof(:, lmp, 1, na, ilSpin)
                        end if
                        this%mat(2, 2 + lo, lp, l, lh) = this%mat(2, 2 + lo, lp, l, lh) &
                                                            & + dot_product(eigVecCoeffs%ccof(m, :, lo, na, ilSpinPr), temp)
                     end do
                  end do

                  do lop = 1, atoms%nlo(itype)
                     lp = atoms%llo(lop, itype)
                     mp = m - mpp
              if ((abs(l - lpp) .le. lp) .and. (lp .le. (l + lpp)) .and. (mod(l + lp + lpp, 2) .eq. 0) .and. (abs(mp) .le. lp)) then
                        fact = cmv*(ImagUnit**(l - lp))*gaunt1(l, lp, lpp, m, mp, mpp, atoms%lmaxd)
                        do na = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype)

                           temp = fact*we(:)*eigVecCoeffs1%ccof(m, :, lo, na, ilSpin)! If not DFPT, this is the base case for rhonmtlo
                           if (l_dfpt) then
                              temp = temp*2.0
                              if (l_gamma) then
                                 temp = temp + fact*we1(:)*eigVecCoeffs%ccof(m, :, lo, na, ilSpin)
                              end if
                           end if
                           this%mat(2 + lo, 2 + lop, l, lp, lh) = this%mat(2 + lo, 2 + lo, l, lp, lh) &
                                     & + dot_product(eigVecCoeffs%ccof(mp, :, lop, na, ilSpinPr), temp)
                        end do

                     end if
                  end do
               end do
            end do
         end do
      end do
   end subroutine rhonmt
end module m_types_denmatrix
