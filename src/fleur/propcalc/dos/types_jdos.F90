!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_jdos
   use m_judft
   use m_types_eigdos
   implicit none
   PRIVATE
   public t_jdos
   TYPE, extends(t_eigdos):: t_jDOS

      REAL, ALLOCATABLE    :: comp(:, :, :, :, :)  !decomposition in percent (neig,0:3,2,n_dos,nkpt)
      REAL, ALLOCATABLE    :: qmtp(:, :, :)      !How much of the state is in the muffin-tin sphere
      REAL, ALLOCATABLE    :: occ(:, :, :)       !Occupation of the j-states
      INTEGER, ALLOCATABLE  :: n_dos_to_na(:)
      ! t2g-based j_eff channels for d-states
      ! comp_jeff_d: (neig, 2, n_dos, nkpt)  index 1=j_eff=1/2, 2=j_eff=3/2 (percent of MT)
      ! comp_jeff_d_mj: (neig, 6, n_dos, nkpt)
      !   indices 1-2: j_eff=1/2, m_j = -1/2,+1/2
      !   indices 3-6: j_eff=3/2, m_j = -3/2,-1/2,+1/2,+3/2
      REAL, ALLOCATABLE    :: comp_jeff_d(:, :, :, :)
      REAL, ALLOCATABLE    :: comp_jeff_d_mj(:, :, :, :)
      REAL, ALLOCATABLE    :: occ_jeff_d(:, :)
      REAL, ALLOCATABLE    :: occ_jeff_d_mj(:, :)

   CONTAINS
      PROCEDURE, PASS :: init => jDOS_init
      PROCEDURE      :: get_weight_eig
      PROCEDURE      :: get_num_weights
      PROCEDURE      :: get_weight_name
      PROCEDURE      :: get_spins
      PROCEDURE      :: calc_jDOS
      procedure      :: postprocessing
   END TYPE t_jDOS
CONTAINS
   subroutine postprocessing(this, noco,nococonv, banddos, alldos, ef)
      use m_types_atoms
      use m_types_noco
      use m_types_nococonv
      use m_types_banddos
      class(t_jDOS), intent(inout):: this
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_banddos), INTENT(IN)    :: banddos
      class(t_eigdos_list), intent(in), optional :: alldos(:)
      real, intent(in), optional :: ef
      return !currently no postprocessing needed for jdos
   end subroutine postprocessing    

   SUBROUTINE calc_jDOS(jDOS, ikpt, noccbd, ev_list, we, atoms, banddos, input, nococonv, itype, radfun, abc_u, abc_d)
      use m_types_atoms
      use m_types_banddos
      use m_types_input
      use m_types_nococonv
      use m_types_radfun
      use m_types_abc
      use m_constants
      use m_clebsch
      CLASS(t_jDOS), INTENT(INOUT)  :: jDOS
      TYPE(t_atoms), INTENT(IN)     :: atoms
      TYPE(t_banddos), INTENT(IN)     :: banddos
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      INTEGER, INTENT(IN)     :: itype
      TYPE(t_radfun), INTENT(IN)     :: radfun
      TYPE(t_abc), INTENT(IN), TARGET :: abc_u, abc_d
      INTEGER, INTENT(IN)     :: ikpt
      INTEGER, INTENT(IN)     :: noccbd
      INTEGER, INTENT(IN)     :: ev_list(:)
      REAL, INTENT(IN)     :: we(:)

      INTEGER, PARAMETER :: lmax = 3 !Maximum l considered in j decomposition
      ! t2g real-harmonic combinations for d-states (l=2):
      !   lm indices in abc%cof: lm = l*(l+1)+m  with l=2 giving lm in [4..8]
      !   m = -2,-1, 0,+1,+2  ->  lm = 4,5,6,7,8
      ! Real spherical harmonics from complex spherical harmonics:
      !   |xy>  = i/sqrt(2)*( Y_2^{-2} + Y_2^{+2} )  -> combine lm4 and lm8
      !            stored as cof difference:
      !              Re part: abc(:,lm4,:,:), Im part: abc(:,lm8,:,:)
      ! Convention: identical to types_orbcomp entries 5-9:
      !   t2g(1) = |xy>  :  i/sqrt(2)*(lm4 + lm8)  -> factor h on lm4,lm8
      !   t2g(2) = |yz>  :  i/sqrt(2)*(lm5 + lm7)  -> orbcomp index 6
      !   t2g(3) = |xz>  :  1/sqrt(2)*(lm5 - lm7)  -> orbcomp index 7
      ! j_eff=3/2 states (m_jt = -3/2,-1/2,+1/2,+3/2) and j_eff=1/2 (m_jt=-1/2,+1/2)
      ! built via CG for l_eff=1 x s=1/2.
      ! Index map into c_mj(1:6): [j_eff=1/2: mj=-1/2,+1/2] [j_eff=3/2: mj=-3/2,-1/2,+1/2,+3/2]

      INTEGER :: n_dos, jcof, icof
      INTEGER :: jjj, iBand, nn, iAtom, iAtom_l, l, jj, j_ind, lmup, lmdown, spin, ilo, ilop
      INTEGER :: imj, i_t2g_l, i_t2g_r, n_r_d
      ! Rotation of reference frame (same feature as in calc_orb_comp)
      TYPE(t_abc), TARGET  :: abc_u_rot, abc_d_rot
      TYPE(t_abc), POINTER :: p_u, p_d
      REAL    :: alpha_use, beta_use, gamma_use
      REAL    :: j,mj, mup, mdown
      REAL    :: facup, facdown, summed, cf
      COMPLEX :: aup, bup, cup, adown, bdown, cdown, cupp, cdownp
      REAL    :: c(0:lmax*2)
      ! t2g j_eff arrays: c_mj(1:6) per band
      REAL    :: c_mj(6), c_jeff(2)
      ! Clebsch-Gordan table for l_eff=1 x s=1/2 -> j_eff:
      ! Rows: (m_leff, m_s) = 6 t2g spinors ordered as
      !   (mt=-1,up),(mt=-1,dn),(mt=0,up),(mt=0,dn),(mt=+1,up),(mt=+1,dn)
      ! mapping to j_eff=1/2 m_j=-1/2,+1/2 and j_eff=3/2 m_j=-3/2,-1/2,+1/2,+3/2
      ! CG(l=1,ml,s=1/2,ms | j,mj):
      !   j=3/2 mj=-3/2: only ml=-1,ms=-1/2  CG=1
      !   j=3/2 mj=-1/2: ml=-1,ms=+1/2 CG=1/sqrt(3); ml=0,ms=-1/2 CG=sqrt(2/3)
      !   j=3/2 mj=+1/2: ml= 0,ms=+1/2 CG=sqrt(2/3); ml=+1,ms=-1/2 CG=1/sqrt(3)
      !   j=3/2 mj=+3/2: only ml=+1,ms=+1/2 CG=1
      !   j=1/2 mj=-1/2: ml=-1,ms=+1/2 CG=-sqrt(2/3); ml=0,ms=-1/2 CG=1/sqrt(3)
      !   j=1/2 mj=+1/2: ml= 0,ms=+1/2 CG=-1/sqrt(3); ml=+1,ms=-1/2 CG=sqrt(2/3)
      REAL, PARAMETER :: r3=sqrt(1.0/3.0), r23=sqrt(2.0/3.0), rh=sqrt(0.5)
      ! t2g spinor basis (6 states):
      !   spinor 1: t2g(|xy>),  spin up
      !   spinor 2: t2g(|xy>),  spin down
      !   spinor 3: t2g(|yz>),  spin up
      !   spinor 4: t2g(|yz>),  spin down
      !   spinor 5: t2g(|xz>),  spin up
      !   spinor 6: t2g(|xz>),  spin down
      ! Map t2g real harmonics to effective m_leff = {-1,0,+1}:
      !   m_leff=-1 <-> |yz>   (spinors 3,4)
      !   m_leff= 0 <-> |xz>   (spinors 5,6)
      !   m_leff=+1 <-> |xy>   (spinors 1,2)
      ! CG coefficients for each (c_mj index, t2g spinor) pair stored as cg_jeff(mj_idx, spinor):
      REAL :: cg_jeff(6, 6)
      ! amplitudes for t2g spinors (ket and bra for 2 radial indices)
      COMPLEX :: amp_ket(6), amp_bra(6)
      COMPLEX :: lc_xy_u, lc_xy_d, lc_yz_u, lc_yz_d, lc_xz_u, lc_xz_d
      COMPLEX :: lc_xy_u2, lc_xy_d2, lc_yz_u2, lc_yz_d2, lc_xz_u2, lc_xz_d2
      INTEGER :: lm_xy_m2, lm_xy_p2, lm_yz_m1, lm_yz_p1, lm_xz_m1, lm_xz_p1
      REAL, PARAMETER :: inv_sqrt2 = sqrt(0.5)

      if (.not. jDOS%l_initialized) return

      ! --- CG table: cg_jeff(i_mj_channel, i_t2g_spinor) ---
      ! channel order: 1=j1/2,mj=-1/2  2=j1/2,mj=+1/2
      !                3=j3/2,mj=-3/2  4=j3/2,mj=-1/2  5=j3/2,mj=+1/2  6=j3/2,mj=+3/2
      ! spinor order: 1=|xy,up> 2=|xy,dn> 3=|yz,up> 4=|yz,dn> 5=|xz,up> 6=|xz,dn>
      ! m_leff mapping: |yz>=m_leff=-1, |xz>=m_leff=0, |xy>=m_leff=+1
      cg_jeff = 0.0
      ! j=1/2, mj=-1/2: CG(l=1,ml=-1,ms=+1/2)=-sqrt(2/3); CG(l=1,ml=0,ms=-1/2)=1/sqrt(3)
      cg_jeff(1, 3) = -r23  ! |yz,up>  ml=-1,ms=+1/2
      cg_jeff(1, 6) =  r3   ! |xz,dn>  ml= 0,ms=-1/2
      ! j=1/2, mj=+1/2: CG(l=1,ml=0,ms=+1/2)=-1/sqrt(3); CG(l=1,ml=+1,ms=-1/2)=sqrt(2/3)
      cg_jeff(2, 5) = -r3   ! |xz,up>  ml= 0,ms=+1/2
      cg_jeff(2, 2) =  r23  ! |xy,dn>  ml=+1,ms=-1/2
      ! j=3/2, mj=-3/2: CG(l=1,ml=-1,ms=-1/2)=1
      cg_jeff(3, 4) =  1.0  ! |yz,dn>  ml=-1,ms=-1/2
      ! j=3/2, mj=-1/2: CG(l=1,ml=-1,ms=+1/2)=1/sqrt(3); CG(l=1,ml=0,ms=-1/2)=sqrt(2/3)
      cg_jeff(4, 3) =  r3   ! |yz,up>  ml=-1,ms=+1/2
      cg_jeff(4, 6) =  r23  ! |xz,dn>  ml= 0,ms=-1/2
      ! j=3/2, mj=+1/2: CG(l=1,ml=0,ms=+1/2)=sqrt(2/3); CG(l=1,ml=+1,ms=-1/2)=1/sqrt(3)
      cg_jeff(5, 5) =  r23  ! |xz,up>  ml= 0,ms=+1/2
      cg_jeff(5, 2) =  r3   ! |xy,dn>  ml=+1,ms=-1/2
      ! j=3/2, mj=+3/2: CG(l=1,ml=+1,ms=+1/2)=1
      cg_jeff(6, 1) =  1.0  ! |xy,up>  ml=+1,ms=+1/2

      ! lm indices for l=2 in abc%cof: lm = l*(l+1)+m
      lm_xy_m2 = 4   ! l=2, m=-2
      lm_xy_p2 = 8   ! l=2, m=+2
      lm_yz_m1 = 5   ! l=2, m=-1
      lm_yz_p1 = 7   ! l=2, m=+1
      lm_xz_m1 = 5   ! l=2, m=-1 (shared with yz; different linear combo)
      lm_xz_p1 = 7   ! l=2, m=+1

      DO iAtom_l = 1, atoms%neq(itype)
         iAtom = iAtom_l - 1 + atoms%firstAtom(itype)
         if (.not. banddos%dos_atom(iAtom)) cycle
         !find index for dos
         DO n_dos = 1, size(banddos%dos_atomlist)
            if (banddos%dos_atomlist(n_dos) == iAtom) exit
         END DO
         ! Rotate abc coefficients into local quantisation frame
         IF (banddos%align_to_spin(iAtom)) THEN
            alpha_use = 0.0
            beta_use  = nococonv%beta(itype)
            IF (ABS(beta_use) > 1.0e-6) THEN
               gamma_use = pi_const/2.0 - nococonv%alph(itype)
            ELSE
               gamma_use = 0.0
            END IF
         ELSE
            alpha_use = banddos%alpha(iAtom)
            beta_use  = banddos%beta(iAtom)
            gamma_use = banddos%gamma(iAtom)
         END IF

         IF (ANY((/alpha_use, beta_use, gamma_use/) .NE. 0.0)) THEN
            abc_u_rot = abc_u%rotate(alpha_use, beta_use, gamma_use, lmax)
            abc_d_rot = abc_d%rotate(alpha_use, beta_use, gamma_use, lmax)
            p_u => abc_u_rot
            p_d => abc_d_rot
         ELSE
            p_u => abc_u
            p_d => abc_d
         END IF
         DO iBand = 1, noccbd
            j_ind = 0
            c = 0.0
            DO l = 0, lmax
            IF (l == 0) THEN
               !s-states (are not split up by SOC)
               DO jjj = 1, radfun%n_r(l)
                  DO jj = 1, radfun%n_r(l)
                   c(0) = c(0) + p_u%cof(iband, 0, jjj, iatom_l)*conjg(p_u%cof(iband, 0, jj, iatom_l))*radfun%integral(jjj, jj, 0, 1, 1)
                   c(0) = c(0) + p_d%cof(iband, 0, jjj, iatom_l)*conjg(p_d%cof(iband, 0, jj, iatom_l))*radfun%integral(jjj, jj, 0, 2, 2)

                  end do
               end do
            ELSE
               DO jj = 1, 2
                  j_ind = j_ind + 1
                  ! j = l +- 1/2
                  j = l + (jj - 1.5)
                  mj = -j
                  DO WHILE (mj <= j)
                     !mj = -l-+1/2, .... , l+-1/2

                     mup = mj - 0.5
                     mdown = mj + 0.5
                     DO icof = 1, radfun%n_r(l)
                        DO jcof = 1, radfun%n_r(l)

                           IF (input%jspins .EQ. 1) THEN
                              mdown = mdown*(-1)
                              spin = 1
                           ELSE
                              spin = 2
                           END IF

                           IF (ABS(mup) <= l) THEN
                              lmup = l*(l + 1) + INT(mup)
                              facup = clebsch(REAL(l), 0.5, mup, 0.5, j, mj)
                              aup = facup*p_u%cof(iBand, lmup, icof, iAtom_l)
                              bup = facup*p_u%cof(iBand, lmup, jcof, iAtom_l)
                           ELSE
                              aup = 0.0
                              bup = 0.0
                           END IF

                           IF (ABS(mdown) <= l) THEN
                              lmdown = l*(l + 1) + INT(mdown)
                              facdown = clebsch(REAL(l), 0.5, mdown, -0.5, j, mj)
                              adown = -1*facdown*p_d%cof(iBand, lmdown, icof, iAtom_l)
                              bdown = -1*facdown*p_d%cof(iBand, lmdown, jcof, iAtom_l)
                           ELSE
                              adown = 0.0
                              bdown = 0.0
                           END IF

                           !c := norm of facup |up> + facdown |down>
                           !We have to write it out explicitely because
                           !of the offdiagonal scalar products that appear
                           c(j_ind) = c(j_ind) &
                                      + aup*CONJG(bup)*radfun%integral(icof, jcof, l, 1, 1) &
                                      + adown*CONJG(bdown)*radfun%integral(icof, jcof, l, 2, 2) &
                                      + aup*CONJG(bdown)*radfun%integral(icof, jcof, l, 1, 2) &
                                      + adown*CONJG(bup)*radfun%integral(icof, jcof, l, 2, 1)
                        end do
                     end do
                     
                     mj = mj + 1
                  END DO
               END DO
            END IF
            END DO
            summed = SUM(c(0:2*lmax))
            IF (summed < 1.0e-30) summed = 1.0e-30
            cf = 100.0/summed
            j_ind = 0
            DO l = 0, 3
            DO jj = 1, 2
               IF (l /= 0) j_ind = j_ind + 1
               jDOS%comp(ev_list(iBand), l, jj, n_dos, ikpt) = c(j_ind)*cf
               jDOS%qmtp(ev_list(iBand), n_dos, ikpt) = 100.0*summed
               jDOS%occ(l, jj, n_dos) = jDOS%occ(l, jj, n_dos) + we(iBand)*c(j_ind)
            END DO
            END DO

            ! --- t2g-based j_eff decomposition for d-states ---
            ! Only meaningful when l=2 radial functions are available
            n_r_d = radfun%n_r(2)
            c_mj = 0.0
            IF (n_r_d > 0) THEN
               DO icof = 1, n_r_d
                  DO jcof = 1, n_r_d
                     ! Build real-harmonic t2g spinor amplitudes from complex abc coefficients.
                     ! Real harmonics expressed in terms of Y_l^m:
                     !   |xy>  =  i/sqrt(2) * (Y_2^{-2} + Y_2^{+2})  [orbcomp entry 5]
                     !   |yz>  =  i/sqrt(2) * (Y_2^{-1} + Y_2^{+1})  [orbcomp entry 6]
                     !   |xz>  =  1/sqrt(2) * (Y_2^{-1} - Y_2^{+1})  [orbcomp entry 7]
                     ! Because abc stores complex Y_l^m coefficients, projecting onto the
                     ! real harmonic requires a phase.  The inner product
                     !   <psi|xy> ~ i/sqrt(2)*(a(lm4)+a(lm8))
                     ! so the amplitude for the t2g orbital part is:
                     lc_xy_u  = CMPLX(0.0,inv_sqrt2)*(p_u%cof(iBand,lm_xy_m2,icof,iAtom_l)+p_u%cof(iBand,lm_xy_p2,icof,iAtom_l))
                     lc_xy_d  = CMPLX(0.0,inv_sqrt2)*(p_d%cof(iBand,lm_xy_m2,icof,iAtom_l)+p_d%cof(iBand,lm_xy_p2,icof,iAtom_l))
                     lc_yz_u  = CMPLX(0.0,inv_sqrt2)*(p_u%cof(iBand,lm_yz_m1,icof,iAtom_l)+p_u%cof(iBand,lm_yz_p1,icof,iAtom_l))
                     lc_yz_d  = CMPLX(0.0,inv_sqrt2)*(p_d%cof(iBand,lm_yz_m1,icof,iAtom_l)+p_d%cof(iBand,lm_yz_p1,icof,iAtom_l))
                     lc_xz_u  = inv_sqrt2*(p_u%cof(iBand,lm_xz_m1,icof,iAtom_l)-p_u%cof(iBand,lm_xz_p1,icof,iAtom_l))
                     lc_xz_d  = inv_sqrt2*(p_d%cof(iBand,lm_xz_m1,icof,iAtom_l)-p_d%cof(iBand,lm_xz_p1,icof,iAtom_l))
                     ! bra (jcof index)
                     lc_xy_u2 = CMPLX(0.0,inv_sqrt2)*(p_u%cof(iBand,lm_xy_m2,jcof,iAtom_l)+p_u%cof(iBand,lm_xy_p2,jcof,iAtom_l))
                     lc_xy_d2 = CMPLX(0.0,inv_sqrt2)*(p_d%cof(iBand,lm_xy_m2,jcof,iAtom_l)+p_d%cof(iBand,lm_xy_p2,jcof,iAtom_l))
                     lc_yz_u2 = CMPLX(0.0,inv_sqrt2)*(p_u%cof(iBand,lm_yz_m1,jcof,iAtom_l)+p_u%cof(iBand,lm_yz_p1,jcof,iAtom_l))
                     lc_yz_d2 = CMPLX(0.0,inv_sqrt2)*(p_d%cof(iBand,lm_yz_m1,jcof,iAtom_l)+p_d%cof(iBand,lm_yz_p1,jcof,iAtom_l))
                     lc_xz_u2 = inv_sqrt2*(p_u%cof(iBand,lm_xz_m1,jcof,iAtom_l)-p_u%cof(iBand,lm_xz_p1,jcof,iAtom_l))
                     lc_xz_d2 = inv_sqrt2*(p_d%cof(iBand,lm_xz_m1,jcof,iAtom_l)-p_d%cof(iBand,lm_xz_p1,jcof,iAtom_l))
                     ! t2g spinor ket/bra arrays (spinor index = orbital x spin):
                     ! 1=|xy,up>  2=|xy,dn>  3=|yz,up>  4=|yz,dn>  5=|xz,up>  6=|xz,dn>
                     amp_ket(1) = lc_xy_u;  amp_ket(2) = lc_xy_d
                     amp_ket(3) = lc_yz_u;  amp_ket(4) = lc_yz_d
                     amp_ket(5) = lc_xz_u;  amp_ket(6) = lc_xz_d
                     amp_bra(1) = lc_xy_u2; amp_bra(2) = lc_xy_d2
                     amp_bra(3) = lc_yz_u2; amp_bra(4) = lc_yz_d2
                     amp_bra(5) = lc_xz_u2; amp_bra(6) = lc_xz_d2

                     ! For each j_eff m_j channel:
                     !   c_mj(ch) += CG_ket . amp_ket  *  conj(CG_bra . amp_bra)  * radial_integral
                     ! Summing over t2g spinors with their CG factors.
                     ! The radial integral depends only on l=2 and spin combination of the two
                     ! spinors; since all t2g spinors share the same l=2 radial function we can
                     ! separate the up-up, dn-dn, up-dn and dn-up pairings.
                     DO imj = 1, 6
                        BLOCK
                           COMPLEX :: ket_u, ket_d, bra_u, bra_d
                           ! project ket/bra onto this m_j channel using CG table
                           ! cg_jeff(imj, s): spinor s has orbital index from the spinor ordering
                           ! Spin-up contribution of the ket: spinors 1,3,5 (odd indices)
                           ket_u = cg_jeff(imj,1)*amp_ket(1) + cg_jeff(imj,3)*amp_ket(3) + cg_jeff(imj,5)*amp_ket(5)
                           ! Spin-dn contribution of the ket: spinors 2,4,6
                           ket_d = cg_jeff(imj,2)*amp_ket(2) + cg_jeff(imj,4)*amp_ket(4) + cg_jeff(imj,6)*amp_ket(6)
                           bra_u = cg_jeff(imj,1)*amp_bra(1) + cg_jeff(imj,3)*amp_bra(3) + cg_jeff(imj,5)*amp_bra(5)
                           bra_d = cg_jeff(imj,2)*amp_bra(2) + cg_jeff(imj,4)*amp_bra(4) + cg_jeff(imj,6)*amp_bra(6)
                           c_mj(imj) = c_mj(imj) &
                                + ket_u*CONJG(bra_u)*radfun%integral(icof,jcof,2,1,1) &
                                + ket_d*CONJG(bra_d)*radfun%integral(icof,jcof,2,2,2) &
                                + ket_u*CONJG(bra_d)*radfun%integral(icof,jcof,2,1,2) &
                                + ket_d*CONJG(bra_u)*radfun%integral(icof,jcof,2,2,1)
                        END BLOCK
                     END DO
                  END DO ! jcof
               END DO ! icof
            END IF ! n_r_d > 0

            ! sum m_j channels to get total j_eff channels
            c_jeff(1) = c_mj(1) + c_mj(2)        ! j_eff=1/2
            c_jeff(2) = c_mj(3) + c_mj(4) + c_mj(5) + c_mj(6) ! j_eff=3/2

            ! normalize by summed MT weight (same normalisation as comp)
            jDOS%comp_jeff_d(ev_list(iBand), 1, n_dos, ikpt) = c_jeff(1)*cf
            jDOS%comp_jeff_d(ev_list(iBand), 2, n_dos, ikpt) = c_jeff(2)*cf
            DO imj = 1, 6
               jDOS%comp_jeff_d_mj(ev_list(iBand), imj, n_dos, ikpt) = c_mj(imj)*cf
            END DO
            ! occupation accumulators
            jDOS%occ_jeff_d(1, n_dos) = jDOS%occ_jeff_d(1, n_dos) + we(iBand)*c_jeff(1)
            jDOS%occ_jeff_d(2, n_dos) = jDOS%occ_jeff_d(2, n_dos) + we(iBand)*c_jeff(2)
            DO imj = 1, 6
               jDOS%occ_jeff_d_mj(imj, n_dos) = jDOS%occ_jeff_d_mj(imj, n_dos) + we(iBand)*c_mj(imj)
            END DO

         END DO ! iBand
      END DO ! iAtom

   END SUBROUTINE calc_jDOS

   pure integer function get_spins(this)
      CLASS(t_jdos), INTENT(IN)::this
      get_spins = 1
   END function

   function get_weight_eig(this, id)
      class(t_jdos), intent(in):: this
      INTEGER, intent(in)      :: id
      real, allocatable:: get_weight_eig(:, :, :)

      integer :: i, l, jj, na, imj

      ALLOCATE (get_weight_eig(size(this%comp, 1), size(this%comp, 5), 1))

      i = 0
      DO na = 1, size(this%comp, 4)
         ! --- legacy 7 channels ---
         DO l = 0, 3
            DO jj = 1, MERGE(1, 2, l == 0)
               i = i + 1
               if (i == id) get_weight_eig(:, :, 1) = this%comp(:, l, jj, na, :)*this%qmtp(:, na, :)/10000.
               if (i > id) RETURN
            END DO
         END DO
         ! --- 2 total j_eff d channels ---
         i = i + 1
         if (i == id) get_weight_eig(:, :, 1) = this%comp_jeff_d(:, 1, na, :)*this%qmtp(:, na, :)/10000.
         if (i > id) RETURN
         i = i + 1
         if (i == id) get_weight_eig(:, :, 1) = this%comp_jeff_d(:, 2, na, :)*this%qmtp(:, na, :)/10000.
         if (i > id) RETURN
         ! --- 6 m_j sub-channels ---
         DO imj = 1, 6
            i = i + 1
            if (i == id) get_weight_eig(:, :, 1) = this%comp_jeff_d_mj(:, imj, na, :)*this%qmtp(:, na, :)/10000.
            if (i > id) RETURN
         END DO
      END DO
   end function

   integer function get_num_weights(this)
      class(t_jdos), intent(in):: this
      ! 7 legacy j-resolved channels + 2 total j_eff + 6 m_j sub-channels = 15 per atom
      get_num_weights = 15*size(this%comp, 4)
   end function

   character(len=20) function get_weight_name(this, id)
      class(t_jdos), intent(in):: this
      INTEGER, intent(in)         :: id
      integer :: i, l, jj, na, imj
      character :: spdfg(0:4) = ["s", "p", "d", "f", "g"]
      character(len=3) :: jname
      ! m_j label table for the 6 d-state j_eff sub-channels
      ! order: j1/2 mj=-1/2,+1/2  then  j3/2 mj=-3/2,-1/2,+1/2,+3/2
      character(len=10), parameter :: mj_names(6) = &
         ["d1-2m-1   ","d1-2m+1   ","d3-2m-3   ","d3-2m-1   ","d3-2m+1   ","d3-2m+3   "]

      i = 0
      DO na = 1, size(this%comp, 4)
         ! --- legacy 7 channels: s, p3/2,p1/2, d3/2,d5/2, f5/2,f7/2 ---
         DO l = 0, 3
            DO jj = -1, MERGE(-1, 1, l == 0), 2
               i = i + 1
               WRITE (jname, '(i1,a,i1)') INT(2*l + jj), '-', 2
               if (i == id) THEN
                  IF (l .EQ. 0) write (get_weight_name, "(a,i0,a)") "jDOS:", this%n_dos_to_na(na), spdfg(l)
                  IF (l .NE. 0) write (get_weight_name, "(a,i0,a,a)") "jDOS:", this%n_dos_to_na(na), spdfg(l), jname
               end if
               if (i > id) RETURN
            END DO
         END DO
         ! --- 2 total j_eff channels for d-states ---
         i = i + 1
         if (i == id) write (get_weight_name, "(a,i0,a)") "jDOS:", this%n_dos_to_na(na), "dj1-2"
         if (i > id) RETURN
         i = i + 1
         if (i == id) write (get_weight_name, "(a,i0,a)") "jDOS:", this%n_dos_to_na(na), "dj3-2"
         if (i > id) RETURN
         ! --- 6 m_j sub-channels for d-states ---
         DO imj = 1, 6
            i = i + 1
            if (i == id) write (get_weight_name, "(a,i0,a)") "jDOS:", this%n_dos_to_na(na), TRIM(mj_names(imj))
            if (i > id) RETURN
         END DO
      END DO

   end function

   SUBROUTINE jDOS_init(thisjDOS, input, banddos, atoms, kpts, eig)

      USE m_types_setup
      USE m_types_kpts

      IMPLICIT NONE

      CLASS(t_jDOS), INTENT(INOUT) :: thisjDOS
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_banddos), INTENT(IN)    :: banddos

      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_kpts), INTENT(IN)    :: kpts
      REAL, INTENT(IN)                      :: eig(:, :, :)
      thisjDOS%l_initialized = .TRUE.
      thisjDOS%n_dos_to_na = banddos%dos_atomlist
      IF (banddos%l_jdos .AND. banddos%dos) THEN
         ALLOCATE (thisjDOS%comp(input%neig, 0:3, 2, size(banddos%dos_atomlist), kpts%nkpt), source=0.0)
         ALLOCATE (thisjDOS%qmtp(input%neig, size(banddos%dos_atomlist), kpts%nkpt), source=0.0)
         ALLOCATE (thisjDOS%occ(0:3, 2, size(banddos%dos_atomlist)), source=0.0)
         ALLOCATE (thisjDOS%comp_jeff_d(input%neig, 2, size(banddos%dos_atomlist), kpts%nkpt), source=0.0)
         ALLOCATE (thisjDOS%comp_jeff_d_mj(input%neig, 6, size(banddos%dos_atomlist), kpts%nkpt), source=0.0)
         ALLOCATE (thisjDOS%occ_jeff_d(2, size(banddos%dos_atomlist)), source=0.0)
         ALLOCATE (thisjDOS%occ_jeff_d_mj(6, size(banddos%dos_atomlist)), source=0.0)
         thisjDOS%eig = eig
      ELSE
         ALLOCATE (thisjDOS%dos(0, 0, 0))
         ALLOCATE (thisjDOS%comp(1, 1, 1, 0, 1), source=0.0)
         ALLOCATE (thisjDOS%qmtp(1, 0, 1), source=0.0)
         ALLOCATE (thisjDOS%occ(1, 1, 0), source=0.0)
         ALLOCATE (thisjDOS%comp_jeff_d(1, 2, 0, 1), source=0.0)
         ALLOCATE (thisjDOS%comp_jeff_d_mj(1, 6, 0, 1), source=0.0)
         ALLOCATE (thisjDOS%occ_jeff_d(2, 0), source=0.0)
         ALLOCATE (thisjDOS%occ_jeff_d_mj(6, 0), source=0.0)
      END IF

      thisjDOS%name_of_dos = "jDOS"

   END SUBROUTINE jDOS_init
end module m_types_jDOS
