!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_desym_potden
   use m_juDFT
   use m_constants
   use m_gaussp
   use m_ylm
   use m_spgrot
   use m_outcdn
   use m_plot, only: getMTSphere
   use m_types_atoms
   use m_types_cell
   use m_types_desym_map
   use m_types_fleurinput
   use m_types_potden
   use m_types_sphhar
   use m_types_stars
   use m_types_sym
   implicit none
   private
   public :: desym_potden, desym_check_potden, potden_on_grid

contains

   subroutine desym_potden(fi, stars, sphhar, fi_nosym, stars_nosym, sphhar_nosym, map, pd, pd_nosym)
      !! Converts a density or potential of the symmetric setup (fi, stars, sphhar) into the desymmetrized setup (fi_nosym, stars_nosym, sphhar_nosym).
      type(t_fleurinput), intent(in)  :: fi, fi_nosym
      type(t_stars),      intent(in)  :: stars, stars_nosym
      type(t_sphhar),     intent(in)  :: sphhar, sphhar_nosym
      type(t_desym_map),  intent(in)  :: map
      type(t_potden),     intent(in)  :: pd
      type(t_potden),     intent(out) :: pd_nosym

      integer :: it_nosym, it

      if (fi%noco%l_noco) call juDFT_error("noco desymmetrization is not implemented yet.", calledby="desym_potden")

      call pd_nosym%init(stars_nosym, fi_nosym%atoms, sphhar_nosym, fi_nosym%vacuum, fi_nosym%noco, fi%input%jspins, pd%potdenType, l_dfpt=allocated(pd%mtIm))
      pd_nosym%iter = pd%iter

      call desym_pw(fi%sym, stars, stars_nosym, pd%pw, pd_nosym%pw)
      if (allocated(pd%pw_w)) then
         allocate(pd_nosym%pw_w, mold=pd_nosym%pw)
         pd_nosym%pw_w = cmplx(0.0, 0.0)
         call desym_pw(fi%sym, stars, stars_nosym, pd%pw_w, pd_nosym%pw_w)
      end if
      if (fi%input%film) call desym_vac(fi, stars, stars_nosym, pd%vac, pd_nosym%vac)

      do it_nosym = 1, fi_nosym%atoms%ntype
         it = map%parent_type(it_nosym)
         call desym_mt(fi, fi_nosym, sphhar, sphhar_nosym, it_nosym, map%type_op(it_nosym), pd%mt(:,:,it,:), pd_nosym%mt(:,:,it_nosym,:))
         if (allocated(pd%mtIm)) call desym_mt(fi, fi_nosym, sphhar, sphhar_nosym, it_nosym, map%type_op(it_nosym), pd%mtIm(:,:,it,:), pd_nosym%mtIm(:,:,it_nosym,:))
         pd_nosym%mtCore(:,it_nosym,:) = pd%mtCore(:,it,:)
         pd_nosym%qint(it_nosym,:)     = pd%qint(it,:)
         pd_nosym%tec(it_nosym,:)      = pd%tec(it,:)
      end do

   end subroutine desym_potden

   subroutine desym_pw(sym, stars, stars_nosym, pw, pw_nosym)
      !! Expands star coefficients of the symmetric setup into the plane-wave coefficients without symmetry (iG: index of G there).
      type(t_sym),   intent(in)    :: sym
      type(t_stars), intent(in)    :: stars, stars_nosym
      complex,       intent(in)    :: pw(:,:)
      complex,       intent(inout) :: pw_nosym(:,:)

      integer :: iStar, iop, iG, nspin
      integer :: gMembers(3,sym%nop)
      complex :: phase, phase_nosym

      nspin = min(size(pw,2), size(pw_nosym,2))
      do iStar = 1, stars%ng3
         call spgrot(sym%nop, sym%symor, sym%mrot, sym%tau, sym%invtab, stars%kv3(:,iStar), gMembers)
         do iop = 1, sym%nop
            iG = stars_nosym%ig(gMembers(1,iop), gMembers(2,iop), gMembers(3,iop))
            if (iG == 0) cycle
            phase       = stars%rgphs(gMembers(1,iop), gMembers(2,iop), gMembers(3,iop))
            phase_nosym = stars_nosym%rgphs(gMembers(1,iop), gMembers(2,iop), gMembers(3,iop))
            if (abs(phase_nosym) < 1e-10) cycle
            pw_nosym(iG,:nspin) = pw(iStar,:nspin) * phase / phase_nosym
         end do
      end do

   end subroutine desym_pw

   subroutine desym_vac(fi, stars, stars_nosym, vac, vac_nosym)
      !! Expands the 2D star coefficients of both vacua into the 2D plane-wave coefficients without symmetry (iG2: index of G_par there).
      type(t_fleurinput), intent(in)    :: fi
      type(t_stars),      intent(in)    :: stars, stars_nosym
      complex,            intent(in)    :: vac(:,:,:,:)
      complex,            intent(inout) :: vac_nosym(:,:,:,:)

      integer :: iStar2, iop, iG2, nspin
      integer :: gMembers(2,fi%sym%nop2)
      complex :: phase, phase_nosym
      complex, allocatable :: vac_both(:,:,:,:)

      ! With nvac = 1 only the upper vacuum is stored; the lower one follows by symmetry
      vac_both = vac
      call stars%fill_2nd_vac(fi%vacuum, vac_both)

      nspin = min(size(vac,4), size(vac_nosym,4))
      do iStar2 = 1, stars%ng2
         call spgrot(fi%sym%nop2, fi%sym%symor, fi%sym%mrot(:2,:2,:), fi%sym%tau(:2,:), fi%sym%invtab, stars%kv2(:,iStar2), gMembers)
         do iop = 1, fi%sym%nop2
            iG2 = stars_nosym%i2g(gMembers(1,iop), gMembers(2,iop))
            if (iG2 == 0) cycle
            phase       = stars%r2gphs(gMembers(1,iop), gMembers(2,iop))
            phase_nosym = stars_nosym%r2gphs(gMembers(1,iop), gMembers(2,iop))
            if (abs(phase_nosym) < 1e-10) cycle
            vac_nosym(:,iG2,:,:nspin) = vac_both(:,iStar2,:,:nspin) * phase / phase_nosym
         end do
      end do

   end subroutine desym_vac

   subroutine desym_mt(fi, fi_nosym, sphhar, sphhar_nosym, it_nosym, iop, mt, mt_nosym)
      !! Rotates the lattice-harmonic expansion of the parent type to the representative atom of type it_nosym
      !! and projects it onto the lattice harmonics without symmetry. The overlap is integrated on a Gauss mesh.
      type(t_fleurinput), intent(in)    :: fi, fi_nosym
      type(t_sphhar),     intent(in)    :: sphhar, sphhar_nosym
      integer,            intent(in)    :: it_nosym, iop
      real,               intent(in)    :: mt(:,0:,:)
      real,               intent(inout) :: mt_nosym(:,0:,:)

      integer :: iatom, nd, nd_nosym, lmax, npts, ip, ilh, ilh_nosym, jri, ispin
      real    :: rot(3,3), x(3)
      real,    allocatable :: vg(:,:), wt(:), lh(:,:), lh_nosym(:,:), ovl(:,:), gram(:,:)
      complex, allocatable :: ylm(:)

      iatom = fi_nosym%atoms%firstAtom(it_nosym)
      nd    = fi%sym%ntypsy(iatom)
      nd_nosym  = fi_nosym%sym%ntypsy(iatom)
      lmax  = max(maxval(sphhar%llh(0:sphhar%nlh(nd),nd)), maxval(sphhar_nosym%llh(0:sphhar_nosym%nlh(nd_nosym),nd_nosym)))
      jri   = fi_nosym%atoms%jri(it_nosym)

      ! Cartesian form of the operation used in outcdn: rho_atom(x) = rho_rep(rot x)
      rot = matmul(fi%cell%amat, matmul(real(fi%sym%mrot(:,:,iop)), fi%cell%bmat)) / tpi_const

      npts = (2*lmax+1) * (lmax+1+mod(lmax+1,2))
      allocate(vg(3,npts), wt(npts), ylm((lmax+1)**2))
      allocate(lh(npts,0:sphhar%nlh(nd)), lh_nosym(npts,0:sphhar_nosym%nlh(nd_nosym)))
      call gaussp(lmax, vg, wt)

      do ip = 1, npts
         call ylm4(lmax, vg(:,ip), ylm)
         call eval_lh(sphhar_nosym, nd_nosym, ylm, lh_nosym(ip,:))
         x = matmul(rot, vg(:,ip))
         call ylm4(lmax, x, ylm)
         call eval_lh(sphhar, nd, ylm, lh(ip,:))
      end do

      allocate(ovl(0:sphhar_nosym%nlh(nd_nosym), 0:sphhar%nlh(nd)), gram(0:sphhar_nosym%nlh(nd_nosym), 0:sphhar_nosym%nlh(nd_nosym)))
      do ilh_nosym = 0, sphhar_nosym%nlh(nd_nosym)
         do ilh = 0, sphhar%nlh(nd)
            ovl(ilh_nosym,ilh) = sum(wt * lh_nosym(:,ilh_nosym) * lh(:,ilh))
         end do
         do ilh = 0, sphhar_nosym%nlh(nd_nosym)
            gram(ilh_nosym,ilh) = sum(wt * lh_nosym(:,ilh_nosym) * lh_nosym(:,ilh))
            if (ilh == ilh_nosym) gram(ilh_nosym,ilh) = gram(ilh_nosym,ilh) - 1.0
         end do
      end do
      if (maxval(abs(gram)) > 1e-8) call juDFT_error("Lattice harmonics of the target group are not orthonormal.", calledby="desym_mt")

      do ispin = 1, min(size(mt,3), size(mt_nosym,3))
         mt_nosym(:jri,:sphhar_nosym%nlh(nd_nosym),ispin) = matmul(mt(:jri,:sphhar%nlh(nd),ispin), transpose(ovl))
      end do

   end subroutine desym_mt

   subroutine eval_lh(sphhar, nd, ylm, lh)
      !! Values of the (real) lattice harmonics of symmetry nd from the spherical harmonics ylm at one point.
      type(t_sphhar), intent(in)  :: sphhar
      integer,        intent(in)  :: nd
      complex,        intent(in)  :: ylm(:)
      real,           intent(out) :: lh(0:)

      integer :: ilh, imem, l, lm
      complex :: s

      do ilh = 0, sphhar%nlh(nd)
         l = sphhar%llh(ilh,nd)
         s = cmplx(0.0, 0.0)
         do imem = 1, sphhar%nmem(ilh,nd)
            lm = l*(l+1) + sphhar%mlh(imem,ilh,nd) + 1
            s = s + sphhar%clnu(imem,ilh,nd) * ylm(lm)
         end do
         lh(ilh) = real(s)
      end do

   end subroutine eval_lh

   subroutine desym_check_potden(fi, stars, sphhar, fi_nosym, stars_nosym, sphhar_nosym, pd, pd_nosym, l_potential, ngrid, maxdev)
      !! Largest deviation between pd and pd_nosym on a regular real-space grid of ngrid^3 points of the unit cell.
      type(t_fleurinput), intent(in)  :: fi, fi_nosym
      type(t_stars),      intent(in)  :: stars, stars_nosym
      type(t_sphhar),     intent(in)  :: sphhar, sphhar_nosym
      type(t_potden),     intent(in)  :: pd, pd_nosym
      logical,            intent(in)  :: l_potential
      integer,            intent(in)  :: ngrid
      real,               intent(out) :: maxdev

      real, allocatable :: vals(:,:), vals_nosym(:,:), points(:,:)

      call potden_on_grid(fi, stars, sphhar, pd, l_potential, ngrid, points, vals)
      call potden_on_grid(fi_nosym, stars_nosym, sphhar_nosym, pd_nosym, l_potential, ngrid, points, vals_nosym)
      maxdev = maxval(abs(vals - vals_nosym))

   end subroutine desym_check_potden

   subroutine potden_on_grid(fi, stars, sphhar, pd, l_potential, ngrid, points, vals)
      !! Values of pd at the ngrid^3 points of a regular grid of the unit cell (cartesian points(3,ngrid^3)), vals(ngrid^3,jspins).
      type(t_fleurinput),   intent(in)  :: fi
      type(t_stars),        intent(in)  :: stars
      type(t_sphhar),       intent(in)  :: sphhar
      type(t_potden),       intent(in)  :: pd
      logical,              intent(in)  :: l_potential
      integer,              intent(in)  :: ngrid
      real, allocatable,    intent(out) :: points(:,:), vals(:,:)

      integer :: ix, iy, iz, ip, jsp, it, iatom, iflag, lv(3)
      real    :: pt(3), p(3)

      allocate(points(3, ngrid**3), vals(ngrid**3, fi%input%jspins))
      ip = 0
      do iz = 0, ngrid-1
         do iy = 0, ngrid-1
            do ix = 0, ngrid-1
               ip = ip + 1
               points(:,ip) = matmul(fi%cell%amat, [real(ix), real(iy), real(iz)] / real(ngrid))
               if (fi%input%film .and. 2*iz >= ngrid) points(:,ip) = points(:,ip) - fi%cell%amat(:,3)
               if (fi%input%film .and. abs(points(3,ip)) > fi%cell%z1) then
                  iflag = 0
                  it = 0
                  iatom = 0
                  lv = 0
                  pt = points(:,ip)
               else
                  ! getMTSphere writes lattvec_index unconditionally, so it must be passed
                  call getMTSphere(fi%input, fi%cell, fi%atoms, points(:,ip), it, iatom, pt, lv)
                  iflag = merge(1, 2, iatom /= 0)
                  if (iatom == 0) pt = points(:,ip)
               end if
               do jsp = 1, fi%input%jspins
                  p = pt
                  call outcdn(p, it, iatom, 1, iflag, jsp, l_potential, stars, fi%vacuum, sphhar, fi%atoms, fi%sym, fi%cell, pd, vals(ip,jsp), lv)
               end do
            end do
         end do
      end do

   end subroutine potden_on_grid

end module m_desym_potden
