!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mt_tofrom_grid
   USE m_types
   PRIVATE
   REAL, PARAMETER    :: d_15 = 1.e-15
   INTEGER, PARAMETER :: ndvgrd = 6 ! this should be consistent across GGA derivative routines
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
      !ENDIF
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
      REAL, ALLOCATABLE :: drm(:,:), drrm(:,:), mm(:,:)
      REAL, ALLOCATABLE :: chlhtot(:,:),chlhdrtot(:,:),chlhdrrtot(:,:)
      INTEGER:: nd, lh, js, jr, kt, k, nsp,j,i,jspV

      !This snippet is crucial to determine over which spins (Only diagonals in colinear case or also off diags in non colin case.)
      IF (noco%l_mtNocoPot) THEN
         jspV=4
      ELSE
         jspV=jspins
      END IF

      nd = sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1)
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
      ENDIF

      !Allocations in mtNoco case
      IF (noco%l_mtNocoPot) THEN
         !General Noco Allocations
         ALLOCATE(mm(atoms%jmtd, 0:sphhar%nlhd))

         !Allocations in case one uses e.g. GGA with mtNoco
         IF (dograds) THEN
            ALLOCATE(drm(atoms%jmtd,0:sphhar%nlhd),drrm(atoms%jmtd, 0:sphhar%nlhd))
            ALLOCATE(chlhtot(atoms%jmtd, 0:sphhar%nlhd),chlhdrtot(atoms%jmtd, 0:sphhar%nlhd),chlhdrrtot(atoms%jmtd, 0:sphhar%nlhd))
         END IF
      END IF

      !Calc magnetization (This is necessary only in mtNoco case)
      IF(noco%l_mtNocoPot) mm(:,:)=SQRT((0.5*(den_mt(:,:,1)-den_mt(:,:,2)))**2+4*den_mt(:,:,3)**2+4*den_mt(:,:,4)**2)


      !Loop to calculate chlh and necessary gradients (if needed)
      DO lh = 0, sphhar%nlh(nd)
         !         calculates gradients of radial charge densities of l=> 0.
         !         rho*ylh/r**2 is charge density. chlh=rho/r**2.
         !         charge density=sum(chlh*ylh).
         !         chlhdr=d(chlh)/dr, chlhdrr=dd(chlh)/drr.

         !Scaling of the magnetic moments in the same way the charge density is scaled in chlh.
         IF(noco%l_mtNocoPot) mm(:,lh)=mm(:,lh)/(atoms%rmsh(:, n)*atoms%rmsh(:, n))

         DO js = 1, jspV
            DO jr = 1, atoms%jri(n)
               chlh(jr, lh, js) = den_mt(jr, lh, js)/(atoms%rmsh(jr, n)*atoms%rmsh(jr, n))
               !Necessary gradients
               IF (dograds) THEN
                  !Colinear case only needs radial derivatives of chlh
                  CALL grdchlh(1, 1, atoms%jri(n), atoms%dx(n), atoms%rmsh(1, n), &
                                                 chlh(1, lh, js), ndvgrd, chlhdr(1, lh, js), chlhdrr(1, lh,js))
                  IF (noco%l_mtNocoPot) THEN
                  !Noco case also needs radial derivatives of mm
                     CALL grdchlh(1, 1, atoms%jri(n), atoms%dx(n), atoms%rmsh(:, n), &
                                                 mm(:,lh), ndvgrd, drm(:,lh), drrm(:,lh))
                 END IF
              END IF
            ENDDO
         ENDDO ! js
      ENDDO   ! lh

      !The following Loop maps chlh on the k-Grid using the lattice harmonics ylh
      kt = 0
      DO jr = 1, atoms%jri(n)
         ! charge density (on extended grid for all jr)
         ! following are at points on jr-th sphere.
         ch_tmp(:, :) = 0.0
         !  generate the densities on an angular mesh (ch_tmp is needed in mkgylm call later on)
         DO js = 1, jspV
            DO lh = 0, sphhar%nlh(nd)
               DO k = 1, nsp
                     ch_tmp(k, js) = ch_tmp(k, js) + ylh(k, lh, nd)*chlh(jr, lh, js)
               ENDDO
            ENDDO
         ENDDO
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

                  !The following snippet maps chlh and its radial derivatives on a colinear system
                  !using mm and its radial derivatives.
                  IF (noco%l_mtNocoPot) THEN
                      IF (js.EQ.1) THEN
                         chlhtot(jr,lh)=0.5*(chlh(jr, lh, 1)+chlh(jr, lh, 2))
                         chlhdrtot(jr,lh)=0.5*(chlhdr(jr, lh, 1)+chlhdr(jr, lh, 2))
                         chlhdrrtot(jr,lh)=0.5*(chlhdrr(jr, lh, 1)+chlhdrr(jr, lh, 2))
                         chlh(jr,lh,js)=chlhtot(jr,lh)+mm(jr,lh)
                         chlhdr(jr, lh, js)=chlhdrtot(jr,lh)+drm(jr,lh)
                         chlhdrr(jr, lh, js)=chlhdrrtot(jr,lh)+drrm(jr,lh)
                      ELSE IF (js.EQ.2) THEN
                         chlh(jr,lh,js)=chlhtot(jr,lh)-mm(jr,lh)
                         chlhdr(jr, lh, js)=chlhdrtot(jr,lh)-drm(jr,lh)
                         chlhdrr(jr, lh, js)=chlhdrrtot(jr,lh)-drrm(jr,lh)
                      ELSE
                         chlh(jr,lh,js)=0
                         chlhdr(jr, lh, js)=0
                         chlhdrr(jr, lh, js)=0
                     END IF
                  END IF

                  !The following loop brings chlhdr and chlhdrr on the k-grid.
                  DO k = 1, nsp
                     chdr(k, js) = chdr(k, js) + ylh(k, lh, nd)*chlhdr(jr, lh, js)
                     chdrr(k, js) = chdrr(k, js) + ylh(k, lh, nd)*chlhdrr(jr, lh, js)
                  ENDDO

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
                  ENDDO
               ENDDO ! lh
            ENDDO   ! js

            !Makegradients 
            IF(jspins>2) CALL mkgylm(2, atoms%rmsh(jr, n), thet, nsp, &
                        ch_tmp, chdr, chdt, chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, grad, kt)
            IF(jspins.LE.2)CALL mkgylm(jspins, atoms%rmsh(jr, n), thet, nsp, &
                        ch_tmp, chdr, chdt, chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, grad, kt)
               ENDIF
         !Set charge to minimum value
         IF (PRESENT(ch)) THEN
            WHERE (ABS(ch_tmp(:nsp,:)) < d_15) ch_tmp(:nsp,:) = d_15
            ch_calc(kt + 1:kt + nsp, :) = ch_tmp(:nsp, :)
         ENDIF
         kt = kt + nsp
      END DO

      !Rotation to local if needed (Indicated by rotch)
      IF (rotch.AND.noco%l_mtNocoPot.AND.(.NOT.dograds)) THEN
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
         ch(:,1:jspins)=ch_calc(:,1:jspins)

      EnD IF
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

      nsp = atoms%nsp()
      nd = sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1)

      DO js = 1, jspins
         !
         kt = 0
         DO jr = 1, atoms%jri(n)
            vpot = v_in(kt + 1:kt + nsp, js)*wt(:)!  multiplicate v_in with the weights of the k-points

            DO lh = 0, sphhar%nlh(nd)
               !
               ! --->        determine the corresponding potential number
               !c            through gauss integration
               !
               vlh = dot_PRODUCT(vpot(:), ylh(:nsp, lh, nd))
               vr(jr, lh, js) = vr(jr, lh, js) + vlh
            ENDDO ! lh
            kt = kt + nsp
         ENDDO   ! jr
      ENDDO

   END SUBROUTINE mt_from_grid

   SUBROUTINE finish_mt_grid()
      DEALLOCATE (ylh, wt, rx, thet, phi)
      IF (ALLOCATED(ylht)) DEALLOCATE (ylht, ylhtt, ylhf, ylhff, ylhtf)
   END SUBROUTINE finish_mt_grid

END MODULE m_mt_tofrom_grid
