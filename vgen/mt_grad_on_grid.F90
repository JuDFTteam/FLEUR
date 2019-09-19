!--------------------------------------------------------------------------------  
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mt_grad_on_grid
   USE m_types
   PRIVATE
   REAL,PARAMETER:: d_15=1.e-15
   INTEGER, PARAMETER :: ndvgrd = 6
   REAL, ALLOCATABLE :: ylh(:, :, :), ylht(:, :, :), ylhtt(:, :, :)
   REAL, ALLOCATABLE :: ylhf(:, :, :), ylhff(:, :, :), ylhtf(:, :, :)
   REAL, ALLOCATABLE :: wt(:), rx(:, :), thet(:), phi(:)
   PUBLIC :: mt_do_grad!, mt_do_div
CONTAINS   
   SUBROUTINE mt_do_grad(jspins, n, atoms, sphhar, sym, den_mt, grad, ch)
      !-----------------------------------------------------------------------------!
      !By use of the radial functions/spherical harmonics representation of a field !
      !component in the Muffin Tins:                                                !
      !                                                                             !
      !Make the gradient of said field in real space and store it.                  !
      !                                                                             !
      !Based on mt_tofrom_grid.F90, A. Neukirchen, September 2019.                  !
      !-----------------------------------------------------------------------------!
      USE m_gaussp
      USE m_lhglptg
      USE m_grdchlh     
      USE m_mkgylm

      IMPLICIT NONE

      INTEGER, INTENT(IN)          :: jspins, n
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)   :: sphhar
      TYPE(t_sym), INTENT(IN)      :: sym
      REAL, INTENT(IN)             :: den_mt(:, 0:, :)
      TYPE(t_gradients), INTENT(INOUT):: grad
      REAL, INTENT(OUT), OPTIONAL  :: ch(:, :)
      
      REAL, ALLOCATABLE :: chlh(:, :, :), chlhdr(:, :, :), chlhdrr(:, :, :)
      REAL, ALLOCATABLE :: ch_tmp(:, :)
      REAL, ALLOCATABLE :: chdr(:, :), chdt(:, :), chdf(:, :) 
      REAL, ALLOCATABLE :: chdrr(:, :), chdtt(:, :), chdff(:, :)
      REAL, ALLOCATABLE :: chdrt(:, :), chdrf(:, :), chdtf(:, :)
      INTEGER:: nd, lh, js, jr, kt, k, nsp
      
      !Generate nspd points on a sherical shell with radius 1.0:
      !The angular mesh is equidistant in phi,
      !theta are the zeros of the legendre polynomials
      ALLOCATE (wt(atoms%nsp()), rx(3, atoms%nsp()), thet(atoms%nsp()), phi(atoms%nsp()))
      CALL gaussp(atoms%lmaxd, rx, wt)
      !Generate the lattice harmonics on the angular mesh:
      ALLOCATE (ylh(atoms%nsp(), 0:sphhar%nlhd, sphhar%ntypsd))
      
      !Generate the derivatives of the lattice harmonics:
      ALLOCATE (ylht, MOLD=ylh)
      ALLOCATE (ylhtt, MOLD=ylh)
      ALLOCATE (ylhf, MOLD=ylh)
      ALLOCATE (ylhff, MOLD=ylh)
      ALLOCATE (ylhtf, MOLD=ylh)
      
      !Calculate the required lattice harmonics and their derivatives:
      CALL lhglptg(sphhar, atoms, rx, atoms%nsp(), xcpot, sym, &
                      ylh, thet, phi, ylht, ylhtt, ylhf, ylhff, ylhtf)
      
      nd = atoms%ntypsy(SUM(atoms%neq(:n - 1)) + 1)
      nsp = atoms%nsp()

      ALLOCATE (ch_tmp(nsp, jspins))
      
      ALLOCATE (chdr(nsp, jspins), chdt(nsp, jspins), chdf(nsp, jspins), &
                chdrr(nsp, jspins),chdtt(nsp, jspins), chdff(nsp, jspins), &
                chdtf(nsp, jspins), chdrt(nsp, jspins), chdrf(nsp, jspins))

      ALLOCATE (chlh(atoms%jmtd, 0:sphhar%nlhd, jspins)) 
      ALLOCATE (chlhdr(atoms%jmtd, 0:sphhar%nlhd, jspins))
      ALLOCATE (chlhdrr(atoms%jmtd, 0:sphhar%nlhd, jspins))
      
      !Calculates the derivatives of the charge density in spherical
      !coordinates as follows [mt are the muffin tin coefficients]: 
      !        
      !chlh=pw/r**2.
      !chlhdr=d(chlh)/dr, chlhdrr=dd(chlh)/drr [a]
      !rho=sum(chlh*ylh)
      ! --> d(rho)/dr=sum(chlhdr*ylh) for radial derivative
      ! --> d(rho)/d[angle]=sum(chlh*ylh[angle]) for the angular ones
      ! --> higher and mixed derivatives follow as one would expect

      !First get chlh, chlhdr and chlhdrr for all grid points [a]:
      DO lh = 0, sphhar%nlh(nd)
         DO js = 1, jspins
            DO jr = 1, atoms%jri(n)
               chlh(jr, lh, js) = den_mt(jr, lh, js)/(atoms%rmsh(jr, n)*atoms%rmsh(jr, n))
            ENDDO ! jr
            CALL grdchlh(1, 1, atoms%jri(n), atoms%dx(n), atoms%rmsh(1, n), chlh(1, lh, js), &
                         ndvgrd, chlhdr(1, lh, js), chlhdrr(1, lh, js))
         ENDDO ! js
      ENDDO ! lh
   
      kt = 0
      DO jr = 1, atoms%jri(n)
         !On extended grid for all jr:
         !Generate rho and all its derivatives on the jr-th sphere:
         ch_tmp(:, :) = 0.0
         chdr(:, :)  = 0.0     ! d(ch)/dr
         chdt(:, :)  = 0.0     ! d(ch)/dtheta
         chdf(:, :)  = 0.0     ! d(ch)/dfai
         chdrr(:, :) = 0.0     ! dd(ch)/drr
         chdtt(:, :) = 0.0     ! dd(ch)/dtt
         chdff(:, :) = 0.0     ! dd(ch)/dff
         chdtf(:, :) = 0.0     ! dd(ch)/dtf
         chdrt(:, :) = 0.0     ! d(d(ch)/dr)dt
         chdrf(:, :) = 0.0     ! d(d(ch)/dr)df
      
         !Calculate [b] and sum [c] the density/derivatives on the angular mesh:
         DO js = 1, jspins
            DO lh = 0, sphhar%nlh(nd)
               DO k = 1, nsp
                  ch_tmp(k, js) = ch_tmp(k, js) + ylh(k, lh, nd)*chlh(jr, lh, js)
                  chdr(k, js) = chdr(k, js) + ylh(k, lh, nd)*chlhdr(jr, lh, js)
                  chdrr(k, js) = chdrr(k, js) + ylh(k, lh, nd)*chlhdrr(jr, lh, js)
                  chdrt(k, js) = chdrt(k, js) + ylht(k, lh, nd)*chlhdr(jr, lh, js)
                  chdrf(k, js) = chdrf(k, js) + ylhf(k, lh, nd)*chlhdr(jr, lh, js)
                  chdt(k, js) = chdt(k, js) + ylht(k, lh, nd)*chlh(jr, lh, js)
                  chdf(k, js) = chdf(k, js) + ylhf(k, lh, nd)*chlh(jr, lh, js)
                  chdtt(k, js) = chdtt(k, js) + ylhtt(k, lh, nd)*chlh(jr, lh, js)
                  chdff(k, js) = chdff(k, js) + ylhff(k, lh, nd)*chlh(jr, lh, js)
                  chdtf(k, js) = chdtf(k, js) + ylhtf(k, lh, nd)*chlh(jr, lh, js)
               ENDDO ! k
            ENDDO ! lh
         ENDDO   ! js
         
         !Combine the derivatives and put them into grad; grad%gr will contain the
         !first derivatives in r, theta, phi on an angular mesh:
         CALL mkgylm(jspins, atoms%rmsh(jr, n), thet, nsp, ch_tmp, chdr, chdt, &
                     chdf, chdrr, chdtt, chdff, chdtf, chdrt, chdrf, grad, kt)
         
         !If requested also pipe out the real space density:
         IF (PRESENT(ch)) THEN
            WHERE (ABS(ch_tmp) < d_15) ch_tmp = d_15
            ch(kt + 1:kt + nsp, :) = ch_tmp(:nsp, :)
         ENDIF
         kt = kt + nsp
      END DO ! jr
   
   END SUBROUTINE mt_do_grad
   
!   SUBROUTINE mt_do_div()
!   Bvec goes in and out as potden
!   mt_do_grad is calles for each component
!   the resulting gradients are combined appropriately into a divergence potden
!   the divergence goes out
! 
!   END SUBROUTINE mt_do_div

END MODULE m_mt_grad_on_grid
