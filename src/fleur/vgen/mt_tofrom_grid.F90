!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mt_tofrom_grid
   USE m_types
   implicit none
   PRIVATE
   REAL, PARAMETER    :: d_15 = 1.e-15
   REAL, ALLOCATABLE :: ylh(:, :, :), ylht(:, :, :), ylhtt(:, :, :)
   REAL, ALLOCATABLE :: ylhf(:, :, :), ylhff(:, :, :), ylhtf(:, :, :)
   REAL, ALLOCATABLE :: wt(:), rx(:, :), thet(:), phi(:)
   PUBLIC :: init_mt_grid, mt_to_grid, mt_from_grid, finish_mt_grid
CONTAINS
   SUBROUTINE init_mt_grid(jspins, atoms, sphhar, dograds, sym, thout, phout)
      USE m_gaussp
      USE m_lhglptg
      USE m_lhglpts
      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: jspins
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)   :: sphhar
      LOGICAL, INTENT(IN)          :: dograds
      TYPE(t_sym), INTENT(IN)      :: sym
      REAL, INTENT(OUT), OPTIONAL  :: thout(:)
      REAL, INTENT(OUT), OPTIONAL  :: phout(:)

      call timestart("init_mt_grid")
      ! generate nspd points on a sherical shell with radius 1.0
      ! angular mesh equidistant in phi,
      ! theta are zeros of the legendre polynomials
      ALLOCATE (wt(atoms%nsp()), rx(3, atoms%nsp()), thet(atoms%nsp()), phi(atoms%nsp()))
      CALL gaussp(atoms%lmaxd, rx, wt)
      ! generate the lattice harmonics on the angular mesh
      ALLOCATE (ylh(atoms%nsp(), 0:sphhar%nlhd, sphhar%ntypsd))
      IF (dograds) THEN
         ALLOCATE (ylht, MOLD=ylh)
         ALLOCATE (ylhtt, MOLD=ylh)
         ALLOCATE (ylhf, MOLD=ylh)
         ALLOCATE (ylhff, MOLD=ylh)
         ALLOCATE (ylhtf, MOLD=ylh)

         CALL lhglptg(sphhar, atoms, rx, atoms%nsp(), dograds, sym, &
                      ylh, thet, phi, ylht, ylhtt, ylhf, ylhff, ylhtf)
         IF (PRESENT(thout)) THEN
            thout=thet
            phout=phi
         END IF

      ELSE
         CALL lhglpts(sphhar, atoms, rx, atoms%nsp(), sym, ylh)
      END IF
      call timestop("init_mt_grid")
   END SUBROUTINE init_mt_grid

   SUBROUTINE mt_to_grid(dograds, jspins, atoms, sym,sphhar,rotch, den_mt, n, noco ,grad, ch)
      USE m_grdchlh
      USE m_mkgylm
      IMPLICIT NONE
      LOGICAL, INTENT(IN)          :: dograds
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_sym), INTENT(IN)      :: sym
      LOGICAL, INTENT(IN)          :: rotch
      TYPE(t_sphhar), INTENT(IN)   :: sphhar
      REAL, INTENT(IN)             :: den_mt(:, 0:, :)
      INTEGER, INTENT(IN)          :: n, jspins
      REAL, INTENT(OUT), OPTIONAL  :: ch(:, :)
      TYPE(t_gradients), INTENT(INOUT):: grad
      TYPE(t_noco), INTENT(IN)     :: noco
      REAL                         :: dentot



      REAL    :: rho_11,rho_22,rho_21r,rho_21i,mx,my,mz,magmom
      REAL    :: rhotot,rho_up,rho_down
      REAL, ALLOCATABLE :: chlh(:, :, :), chlhdr(:, :, :), chlhdrr(:, :, :)
      REAL, ALLOCATABLE :: chdr(:, :), chdt(:, :), chdf(:, :), ch_tmp(:, :),ch_calc(:,:)
      REAL, ALLOCATABLE :: chdrr(:, :), chdtt(:, :), chdff(:, :), chdtf(:, :)
      REAL, ALLOCATABLE :: chdrt(:, :), chdrf(:, :)
      REAL, ALLOCATABLE :: ch_ud(:, :)
      INTEGER:: nd, lh, js, jr, kt, k, nsp,j,i,jspV
      LOGICAL :: l_nocoden

      call timestart("mt_to_grid")
      !Full 2x2 density matrix only if the field carries all four components
      l_nocoden = any(noco%l_unrestrictMT) .AND. SIZE(den_mt,3)>=4
      jspV = MERGE(4, jspins, l_nocoden)

      nd = sym%ntypsy(atoms%firstAtom(n))
      nsp = atoms%nsp()

      !General Allocations
      ALLOCATE (chlh(atoms%jmtd, 0:sphhar%nlhd, jspV))
      ALLOCATE (ch_tmp(nsp, jspV),ch_calc(nsp*atoms%jmtd, jspV))

      !Allocations in dograds case
      IF (dograds) THEN
         ALLOCATE (chdr(nsp, jspV), chdt(nsp, jspV), chdf(nsp, jspV), chdrr(nsp, jspV), &
                   chdtt(nsp, jspV), chdff(nsp, jspV), chdtf(nsp, jspV), chdrt(nsp, jspV), &
                   chdrf(nsp, jspV))
         ALLOCATE (chlhdr(atoms%jmtd, 0:sphhar%nlhd, jspV))
         ALLOCATE (chlhdrr(atoms%jmtd, 0:sphhar%nlhd, jspV))
         IF (l_nocoden) ALLOCATE (ch_ud(nsp, 2))
      END IF

      !Loop to calculate chlh and necessary gradients (if needed)
      DO lh = 0, sphhar%nlh(nd)
         !         calculates gradients of radial charge densities of l=> 0.
         !         rho*ylh/r**2 is charge density. chlh=rho/r**2.
         !         charge density=sum(chlh*ylh).
         !         chlhdr=d(chlh)/dr, chlhdrr=dd(chlh)/drr.

         DO js = 1, jspV
            chlh(1:atoms%jri(n), lh, js) = den_mt(1:atoms%jri(n), lh, js) / &
                                           (atoms%rmsh(1:atoms%jri(n), n)*atoms%rmsh(1:atoms%jri(n), n))

            !Necessary gradients
            IF (dograds) THEN
               !Colinear case only needs radial derivatives of chlh
               CALL grdchlh(atoms%dx(n), chlh(1:atoms%jri(n), lh, js),  chlhdr(1:, lh, js), chlhdrr(1:, lh,js),atoms%rmsh(:,n))
            END IF
         END DO ! js
      END DO   ! lh

      !The following Loop maps chlh on the k-Grid using the lattice harmonics ylh
      !$OMP parallel do default( NONE ) &
      !$OMP SHARED(atoms,sphhar,noco,n,nsp,jspV,nd,ylh,chlh,dograds,chlhdr,chlhdrr,l_nocoden)&
      !$OMP SHARED(ylht,ylhf,ylhtt,ylhff,ylhtf,jspins,thet,grad,ch,ch_calc) &
      !$OMP private(kt,ch_tmp,js,lh,k,chdr,chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf,ch_ud)
      DO jr = 1, atoms%jri(n)
         kt = (jr-1)*nsp
         ! charge density (on extended grid for all jr)
         ! following are at points on jr-th sphere.
         ch_tmp(:, :) = 0.0
         !  generate the densities on an angular mesh (ch_tmp is needed in mkgylm call later on)
         DO js = 1, jspV
            DO lh = 0, sphhar%nlh(nd)
               DO k = 1, nsp
                     ch_tmp(k, js) = ch_tmp(k, js) + ylh(k, lh, nd)*chlh(jr, lh, js)
               END DO
            END DO
         END DO

         !Initialize derivatives of ch on grid if needed.
         IF (dograds) THEN
            chdr(:, :) = 0.0     ! d(ch)/dr
            chdt(:, :) = 0.0     ! d(ch)/dtheta
            chdf(:, :) = 0.0     ! d(ch)/dfai
            chdrr(:, :) = 0.0     ! dd(ch)/drr
            chdtt(:, :) = 0.0     ! dd(ch)/dtt
            chdff(:, :) = 0.0     ! dd(ch)/dff
            chdtf(:, :) = 0.0     ! dd(ch)/dtf
            chdrt(:, :) = 0.0     ! d(d(ch)/dr)dt
            chdrf(:, :) = 0.0     ! d(d(ch)/dr)df
            !  generate the derivatives on an angular mesh
            DO js = 1, jspV
               DO lh = 0, sphhar%nlh(nd)

                  !The following loop brings chlhdr and chlhdrr on the k-grid.
                  DO k = 1, nsp
                     chdr(k, js) = chdr(k, js) + ylh(k, lh, nd)*chlhdr(jr, lh, js)
                     chdrr(k, js) = chdrr(k, js) + ylh(k, lh, nd)*chlhdrr(jr, lh, js)
                  END DO

                  !This loop calculates the other derviatives of ch (Angular terms) on the k-grid
                  !by using the lattice harmonics derivatives and chlh with its derivatives.
                  DO k = 1, nsp
                     chdrt(k, js) = chdrt(k, js) + ylht(k, lh, nd)*chlhdr(jr, lh, js)
                     chdrf(k, js) = chdrf(k, js) + ylhf(k, lh, nd)*chlhdr(jr, lh, js)
                     chdt(k, js) = chdt(k, js) + ylht(k, lh, nd)*chlh(jr, lh, js)
                     chdf(k, js) = chdf(k, js) + ylhf(k, lh, nd)*chlh(jr, lh, js)
                     chdtt(k, js) = chdtt(k, js) + ylhtt(k, lh, nd)*chlh(jr, lh, js)
                     chdff(k, js) = chdff(k, js) + ylhff(k, lh, nd)*chlh(jr, lh, js)
                     chdtf(k, js) = chdtf(k, js) + ylhtf(k, lh, nd)*chlh(jr, lh, js)
                  END DO
               END DO ! lh
            END DO   ! js
         !Rotation to local if needed (Indicated by rotch)
            !Makegradients
            !IF(jspins>2) CALL mkgylm(2, atoms%rmsh(jr, n), thet, nsp, &
            !            ch_tmp, chdr, chdt, chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, grad, kt)
            !IF(jspins.LE.2)
            IF (l_nocoden) THEN
               !Gradients of rho_up/down=(rho+-|m|)/2 in the local frame
               CALL noco_to_updown(ch_tmp, chdr, chdt, chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, ch_ud)
               CALL mkgylm(2, atoms%rmsh(jr, n), thet, nsp, &
                           ch_ud, chdr, chdt, chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, grad, kt)
            ELSE
               CALL mkgylm(jspins, atoms%rmsh(jr, n), thet, nsp, &
                           ch_tmp, chdr, chdt, chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, grad, kt)
            END IF
         END IF
         !Set charge to minimum value
         IF (PRESENT(ch)) THEN
            WHERE (ABS(ch_tmp(:nsp,:)) < d_15) ch_tmp(:nsp,:) = d_15
            ch_calc(kt + 1:kt + nsp, :) = ch_tmp(:nsp, :)
         END IF
      END DO
      !$OMP END PARALLEL DO


      IF (PRESENT(ch)) THEN
      !Rotation to local if needed (Indicated by rotch)
         IF (rotch.AND.l_nocoden) THEN
            DO jr = 1,nsp*atoms%jri(n)
               rho_11  = ch_calc(jr,1)
               rho_22  = ch_calc(jr,2)
               rho_21r = ch_calc(jr,3)
               rho_21i = ch_calc(jr,4)
               mx      =  2*rho_21r
               my      = -2*rho_21i
               mz      = (rho_11-rho_22)
               magmom  = SQRT(mx**2 + my**2 + mz**2)
               rhotot  = rho_11 + rho_22
               rho_up  = (rhotot + magmom)/2
               rho_down= (rhotot - magmom)/2
               ch(jr,1) = rho_up
               ch(jr,2) = rho_down
            END DO
         ELSE
            ch(:nsp*atoms%jri(n),1:jspins)=ch_calc(:nsp*atoms%jri(n),1:jspins)

         END IF
      END IF
      call timestop("mt_to_grid")
   END SUBROUTINE mt_to_grid

   SUBROUTINE mt_from_grid(atoms, sym, sphhar, n, jspins, v_in, vr)
      IMPLICIT NONE
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_sphhar), INTENT(IN):: sphhar
      INTEGER, INTENT(IN)       :: jspins, n
      REAL, INTENT(IN)          :: v_in(:, :)
      REAL, INTENT(INOUT)       :: vr(:, 0:, :)

      REAL    :: vpot(atoms%nsp()), vlh
      INTEGER :: js, kt, lh, jr, nd, nsp

      call timestart("mt_from_grid")

      nsp = atoms%nsp()
      nd = sym%ntypsy(atoms%firstAtom(n))

      DO js = 1, jspins
         !$OMP parallel do default( NONE ) &
         !$OMP SHARED(atoms,sphhar,n,v_in,nsp,wt,nd,ylh,vr,js)&
         !$OMP private(vpot,kt,lh,vlh)
         DO jr = 1, atoms%jri(n)
            kt = (jr-1)*nsp
            vpot = v_in(kt + 1:kt + nsp, js)*wt(:)!  multiplicate v_in with the weights of the k-points

            DO lh = 0, sphhar%nlh(nd)
               !
               ! --->        determine the corresponding potential number
               !c            through gauss integration
               !
               vlh = dot_PRODUCT(vpot(:), ylh(:nsp, lh, nd))
               vr(jr, lh, js) = vr(jr, lh, js) + vlh
            END DO ! lh
         END DO   ! jr
         !$OMP END PARALLEL DO
      END DO
      call timestop("mt_from_grid")
   END SUBROUTINE mt_from_grid

   SUBROUTINE noco_to_updown(ch, chdr, chdt, chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, ch_ud)
      !Converts the 2x2 density matrix (11,22,Re21,Im21) and its derivatives on the grid
      !into rho_up/down=(rho+-|m|)/2. The derivatives are overwritten in columns 1:2.
      IMPLICIT NONE
      REAL, INTENT(IN)    :: ch(:, :)
      REAL, INTENT(INOUT) :: chdr(:, :), chdt(:, :), chdf(:, :)
      REAL, INTENT(INOUT) :: chdrr(:, :), chdtt(:, :), chdff(:, :), chdtf(:, :), chdrt(:, :), chdrf(:, :)
      REAL, INTENT(OUT)   :: ch_ud(:, :)

      REAL, PARAMETER    :: eps = 1.e-10
      !second derivatives rr,tt,ff,tf,rt,rf as pairs of first derivatives r,t,f
      INTEGER, PARAMETER :: ia(6) = [1, 2, 3, 2, 1, 1], ib(6) = [1, 2, 3, 3, 2, 3]
      REAL    :: m(3), d1(3, 3), d2(3, 6), am, dam(3), ddam(6)
      INTEGER :: k, i

      DO k = 1, SIZE(ch, 1)
         m = mvec(ch(k, :))
         d1(:, 1) = mvec(chdr(k, :)); d1(:, 2) = mvec(chdt(k, :)); d1(:, 3) = mvec(chdf(k, :))
         d2(:, 1) = mvec(chdrr(k, :)); d2(:, 2) = mvec(chdtt(k, :)); d2(:, 3) = mvec(chdff(k, :))
         d2(:, 4) = mvec(chdtf(k, :)); d2(:, 5) = mvec(chdrt(k, :)); d2(:, 6) = mvec(chdrf(k, :))
         am = NORM2(m)
         IF (am > eps) THEN
            DO i = 1, 3
               dam(i) = DOT_PRODUCT(m, d1(:, i))/am
            END DO
            DO i = 1, 6
               ddam(i) = (DOT_PRODUCT(d1(:, ia(i)), d1(:, ib(i))) + DOT_PRODUCT(m, d2(:, i)))/am &
                         - dam(ia(i))*dam(ib(i))/am
            END DO
         ELSE
            dam = 0.0; ddam = 0.0
         END IF
         ch_ud(k, 1) = 0.5*(ch(k, 1) + ch(k, 2) + am)
         ch_ud(k, 2) = 0.5*(ch(k, 1) + ch(k, 2) - am)
         CALL updown(chdr(k, :), dam(1));   CALL updown(chdt(k, :), dam(2));   CALL updown(chdf(k, :), dam(3))
         CALL updown(chdrr(k, :), ddam(1)); CALL updown(chdtt(k, :), ddam(2)); CALL updown(chdff(k, :), ddam(3))
         CALL updown(chdtf(k, :), ddam(4)); CALL updown(chdrt(k, :), ddam(5)); CALL updown(chdrf(k, :), ddam(6))
      END DO
   CONTAINS
      PURE FUNCTION mvec(c)
         REAL, INTENT(IN) :: c(:)
         REAL             :: mvec(3)
         mvec = [2*c(3), -2*c(4), c(1) - c(2)]
      END FUNCTION mvec

      PURE SUBROUTINE updown(c, dm)
         REAL, INTENT(INOUT) :: c(:)
         REAL, INTENT(IN)    :: dm
         REAL                :: drho
         drho = c(1) + c(2)
         c(1) = 0.5*(drho + dm)
         c(2) = 0.5*(drho - dm)
      END SUBROUTINE updown
   END SUBROUTINE noco_to_updown

   SUBROUTINE finish_mt_grid()
      implicit NONE
      call timestart("finish_mt_grid")
      DEALLOCATE (ylh, wt, rx, thet, phi)
      IF (ALLOCATED(ylht)) DEALLOCATE (ylht, ylhtt, ylhf, ylhff, ylhtf)
      call timestop("finish_mt_grid")
   END SUBROUTINE finish_mt_grid

END MODULE m_mt_tofrom_grid
