! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_cdnmt
   !! Archived comment:
   !!
   !! This subroutine calculates the spherical and non-spherical charge-
   !! density and the orbital moment inside the muffin-tin spheres.
   !! Philipp Kurz 2000-02-03
contains

   subroutine cdnmt(jspd, input, atoms, sym, sphhar, noco, jsp_start, jsp_end, enpara, banddos, &
                    vr, denCoeffs, usdus, orb, denCoeffsOffdiag, rho, hub1inp, moments, jDOS, hub1data, rhoIm)
      !! Current situation:
      !!
      !! This routine calculates density contributions
      !! $$\rho_{L}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha}(r)=
      !! \sum_{l',l,\lambda',\lambda,s}d_{l',l,L,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha}
      !! u_{l',\lambda',s}^{\sigma_{\alpha}',\alpha}(r)u_{l,\lambda,s}^{\sigma_{\alpha},\alpha}(r)$$
      !! \(s\) is the index for the big/small components yielded by the
      !! scalar-relativistic Schrödinger equation.

      use m_types
      use m_constants
      use m_cdnmtlo
      use m_radfun
      use m_orbmom2
      use m_xmlOutput
      use m_types_orbcomp
      use m_types_jDOS
      use m_types_mcd
      use m_intgr

      implicit none

      type(t_input), intent(IN)    :: input
      type(t_usdus), intent(INOUT) :: usdus !in fact only the lo part is intent(in)
      type(t_noco), intent(IN)    :: noco
      type(t_sphhar), intent(IN)    :: sphhar
      type(t_atoms), intent(IN)    :: atoms
      type(t_sym), intent(IN)    :: sym
      type(t_enpara), intent(IN)    :: enpara
      type(t_banddos), intent(IN)    :: banddos
      type(t_hub1inp), intent(IN)    :: hub1inp
      type(t_orb), intent(IN)    :: orb
      type(t_denCoeffs), intent(IN)    :: denCoeffs
      type(t_denCoeffsOffdiag), intent(IN)    :: denCoeffsOffdiag

      type(t_jDOS), optional, intent(IN)    :: jDOS
      type(t_moments), optional, intent(INOUT) :: moments
      type(t_hub1data), optional, intent(INOUT) :: hub1data

      integer, intent(IN) :: jsp_start, jsp_end, jspd

      real, intent(IN) :: vr(atoms%jmtd, atoms%ntype, jspd)
      real, intent(INOUT) :: rho(:, 0:, :, :)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,jspd)

      real, optional, intent(INOUT) :: rhoIm(:, 0:, :, :)

      integer, parameter :: lcf = 3

      integer :: itype, na, nd, l, m, lp, llp, lh, j, ispin, noded, nodeu, llpb, natom, jj, n_dos
      integer :: ilo, ilop, i, i_hia, i_exc, icoef, jcoef
      real    :: wronk, qmtt, radThomson, tempFactor, sumlm
      complex :: cs, rho21
      logical :: l_hia, l_performSpinavg

      real              :: qmtl(0:atoms%lmaxd, jspd, atoms%ntype), qmtllo(0:atoms%lmaxd), vrTmp(atoms%jmtd)
      real, allocatable  :: vr0(:, :)
      character(LEN=20) :: attributes(6)
      real              :: denThomson(atoms%jmtd, -1:3)
      real              :: overlapR3(atoms%jmtd, 2, 2)
      real              :: ovlpInt(2, 2)

      real, allocatable :: Rf(:, :, :, :, :)

      allocate (Rf(atoms%jmtd, 2, 0:1, 0:atoms%lmaxd, minval(ispin, ispinpr):maxval(ispin, ispinpr)))

      vrTmp = vr0(:, ispin)!check if spin averaged potential should be used
      call radfun(l, itype, ispin, enpara%el0(l, itype, ispin), vrTmp, atoms, &
                  Rf(1, 1, 0, l, ispin), Rf(1, 1, 1, l, ispin), usdus, nodeu, noded, wronk)
      !construct second radial function if needed
      if (ispin /= ispinpr) call radfun(l, itype, ispinpr, enpara%el0(l, itype, ispinpr), vrTmp, atoms, &
                                        Rf(1, 1, 0, l, ispinpr), Rf(1, 1, 1, l, ispinpr), usdus, nodeu, noded, wronk)

      !Now LO part is needed
      call timestart("cdnmt LO")
      call cdnmtlo(itype, ispinpr, ispin, input, atoms, sphhar, sym, usdus, orb, noco, &
                   enpara%ello0(:, itype, :), vr0(:, :), denCoeffs, &
                   Rf(:, :, 0, 0:, ispinpr), Rf(:, :, 1, 0:, ispinpr), &
                   rho(:, 0:), moments=moments, &
                   rhoIm=rhoIm(:, 0:), f2=Rf(:, :, 0, 0:, ispin), g2=Rf(:, :, 1, 0:, ispin))
      call timestop("cdnmt LO")

      do lh = 0, sphhar%nlh(sym%ntypsy(atoms%firstAtom(itype)))
         do l = 0, atoms%lmax(itype)
            do lp = 0, merge(l, atoms%lmax(itype), present(moments))
               llp = (l*(l + 1))/2 + lp
               if (.not. present(moments)) llp = lp*(atoms%lmax(itype) + 1) + l
               if (lh > 0 .and. atoms%l_outputCFpot(itype) .and. atoms%l_outputCFremove4f(itype) &
                   .and. (l == lcf .and. lp == lcf)) cycle !Exclude non-spherical contributions for CF

               do j = 1, atoms%jri(itype)
                  cs = 0.0
                  do icoef = lbound(denCoeffs%nmt_coeff, 4), ubound(denCoeffs%nmt_coeff, 4) !Loop over radial functions (usually 0 and 1)
                     do jcoef = lbound(denCoeffs%nmt_coeff, 5), ubound(denCoeffs%nmt_coeff, 5) !Loop over radial functions (usually 0 and 1)
                            cs = cs+ denCoeffs%nmt_coeff(llp,lh,itype,icoef,jcoef,ispinpr,ispin)*(Rf(j,1,icoef,lp,ispin)*Rf(j,1,jcoef,l,ispinpr)+ Rf(j,2,icoef,lp,ispin)*Rf(j,2,jcoef,l,ispinpr))
                     end do
                  end do
                  rho(j, lh) = rho(j, lh) + real(cs)/atoms%neq(itype)
                  if ((l <= input%lResMax) .and. (lp <= input%lResMax) .and. ispin == ispinpr .and. present(moments)) &
                   moments%rhoLRes(j, lh, llp, itype, ispin) = moments%rhoLRes(j, lh, llp, itype, ispin) + real(cs)/atoms%neq(itype)
                  if (present(rhoIm)) rhoIm(j, lh) = rhoIm(j, lh) + aimag(cs)/atoms%neq(itype)
                  end if
               end do
            end do
         end do
         end do

         end subroutine cdnmt
      end module m_cdnmt
