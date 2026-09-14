!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vis_xc
   USE m_juDFT
   use m_convol
   !     ******************************************************
   !     subroutine generates the exchange-correlation potential
   !     in the interstitial region    c.l.fu
   !     including gradient corrections. t.a. 1996.
   !     ******************************************************
CONTAINS
   SUBROUTINE dfpt_vis_xc(stars,starsq,sym,cell,den,den1,xcpot,input,vTot)

      USE m_pw_tofrom_grid
      USE m_types
      USE m_constants
      USE m_types_xcpot_libxc
      USE m_dfpt_gga_kernel
      IMPLICIT NONE

      CLASS(t_xcpot),INTENT(IN)     :: xcpot
      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_stars),INTENT(IN)      :: stars, starsq
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_potden),INTENT(IN)     :: den, den1
      TYPE(t_potden),INTENT(INOUT)  :: vTot

      TYPE(t_gradients) :: gradRho, gradRho1, gradRho1Im
      TYPE(t_gradients) :: gradDrivsigma, gradDrivsigma1, gradDrivsigma1Im
      TYPE(t_potden) :: vTotIm

      REAL, ALLOCATABLE :: rho(:,:), rho1(:,:), rho1Im(:,:)
      REAL, ALLOCATABLE :: v_xc1(:,:),v_xc1Im(:,:),f_xc(:,:)
      !Kernels come out of libxc with the spin-like index first; drivsigmaT is the
      !transpose needed to put drivsigma back on the grid.
      REAL, ALLOCATABLE :: drivsigma(:,:),driv2rho2(:,:),driv2rhosigma(:,:),driv2sigma2(:,:)
      REAL, ALLOCATABLE :: drivsigmaT(:,:),drivsigma1(:,:),drivsigma1Im(:,:)
      COMPLEX, ALLOCATABLE :: drivsigmaPW(:,:),drivsigma1PW(:,:),drivsigma1PWIm(:,:)
      LOGICAL, ALLOCATABLE :: l_cut(:)
      INTEGER           :: iSpin, jSpin, fxcSpin, i, nfxc, n_sigma, npt

      REAL, PARAMETER   :: den_cut = 1.E-6

      nfxc = 2 * input%jspins - 1
      n_sigma = MERGE(1,3,input%jspins==1)

      IF (ALLOCATED(vTotIm%pw)) DEALLOCATE(vTotIm%pw)
      ALLOCATE(vTotIm%pw,mold=vTot%pw)
      vTotIm%pw = CMPLX(0.0,0.0)

      call timestart("init_pw_grid")
      CALL init_pw_grid(stars,sym,cell,xcpot)
      call timestop("init_pw_grid")

      !Put the charge on the grid, in GGA case also calculate gradients.
      !den1 lives on the q-shifted stars, so its gradient carries i(q+G) and is complex.
      call timestart("pw_to_grid")
      CALL pw_to_grid(xcpot%needs_grad(),input%jspins,.FALSE.,stars,cell,den%pw,gradRho,xcpot,rho)
      IF (xcpot%needs_grad()) THEN
         CALL pw_to_grid(.TRUE.,input%jspins,.FALSE.,starsq,cell,den1%pw,gradRho1,xcpot,rho1,rho1Im,gradim=gradRho1Im)
      ELSE
         CALL pw_to_grid(.FALSE.,input%jspins,.FALSE.,starsq,cell,den1%pw,gradRho1,xcpot,rho1,rho1Im)
      END IF
      call timestop("pw_to_grid")

      npt = SIZE(rho,1)
      ALLOCATE(v_xc1,mold=rho)
      ALLOCATE(v_xc1Im,mold=rho)

      IF (xcpot%needs_grad()) THEN
         ALLOCATE(l_cut(npt))
         l_cut = ANY(ABS(rho)<den_cut,dim=2)

         call timestart("apply_cutoffs")
         CALL xcpot%apply_cutoffs(den_cut,rho,gradRho)
         call timestop("apply_cutoffs")

         ALLOCATE(drivsigma(n_sigma,npt),driv2rho2(nfxc,npt))
         ALLOCATE(driv2rhosigma(input%jspins*n_sigma,npt),driv2sigma2(MERGE(1,6,input%jspins==1),npt))
         SELECT TYPE(xcpot)
         TYPE IS (t_xcpot_libxc)
            CALL xcpot%get_fxc_gga(input%jspins,rho,gradRho%sigma,drivsigma,driv2rho2,driv2rhosigma,driv2sigma2)
         END SELECT

         ALLOCATE(drivsigma1(npt,n_sigma),drivsigma1Im(npt,n_sigma))
         call timestart("dfpt_gga_assemble")
         CALL dfpt_gga_assemble(drivsigma,driv2rho2,driv2rhosigma,driv2sigma2,gradRho,gradRho1,rho1,v_xc1,drivsigma1)
         CALL dfpt_gga_assemble(drivsigma,driv2rho2,driv2rhosigma,driv2sigma2,gradRho,gradRho1Im,rho1Im,v_xc1Im,drivsigma1Im)
         call timestop("dfpt_gga_assemble")

         !drivsigma1 is expanded in stars below, so a spike in the tail of the
         !density would pollute the whole cell.
         DO i = 1, npt
            IF (l_cut(i)) THEN
               drivsigma1(i,:) = 0.0
               drivsigma1Im(i,:) = 0.0
            END IF
         END DO

         call timestart("gradient of drivsigma")
         !Ground state drivsigma; unshifted stars and real.
         ALLOCATE(drivsigmaT(npt,n_sigma),drivsigmaPW(stars%ng3,n_sigma))
         drivsigmaT = TRANSPOSE(drivsigma)
         drivsigmaPW = CMPLX(0.0,0.0)
         CALL pw_from_grid(stars,drivsigmaT,drivsigmaPW)
         CALL pw_to_grid(.TRUE.,n_sigma,.FALSE.,stars,cell,drivsigmaPW,gradDrivsigma,xcpot)

         !Response of drivsigma; q-shifted stars, so both channels are needed.
         ALLOCATE(drivsigma1PW(starsq%ng3,n_sigma),drivsigma1PWIm(starsq%ng3,n_sigma))
         drivsigma1PW = CMPLX(0.0,0.0)
         drivsigma1PWIm = CMPLX(0.0,0.0)
         CALL pw_from_grid(starsq,drivsigma1,drivsigma1PW)
         CALL pw_from_grid(starsq,drivsigma1Im,drivsigma1PWIm)
         drivsigma1PW = drivsigma1PW + ImagUnit * drivsigma1PWIm
         CALL pw_to_grid(.TRUE.,n_sigma,.FALSE.,starsq,cell,drivsigma1PW,gradDrivsigma1,xcpot,gradim=gradDrivsigma1Im)
         call timestop("gradient of drivsigma")

         call timestart("dfpt_gga_grdotgr")
         CALL dfpt_gga_grdotgr(gradRho,gradRho1,gradDrivsigma,gradDrivsigma1,v_xc1)
         CALL dfpt_gga_grdotgr(gradRho,gradRho1Im,gradDrivsigma,gradDrivsigma1Im,v_xc1Im)
         call timestop("dfpt_gga_grdotgr")
      ELSE
         ALLOCATE(f_xc(npt,nfxc))
         CALL xcpot%get_fxc_lda(input%jspins, rho, f_xc)

         v_xc1 = 0.0
         v_xc1Im = 0.0
         DO iSpin = 1, input%jspins
             DO jSpin = 1, input%jspins
                 fxcSpin = iSpin + jSpin - 1
                 v_xc1(:, iSpin) = v_xc1(:, iSpin) + f_xc(:, fxcSpin) * rho1(:, jSpin)
                 v_xc1Im(:, iSpin) = v_xc1Im(:, iSpin) + f_xc(:, fxcSpin) * rho1Im(:, jSpin)
             END DO
         END DO
      END IF

      !Put the potentials in rez. space.
      call timestart("pw_from_grid")
      CALL  pw_from_grid(starsq,v_xc1,vTot%pw)
      CALL  pw_from_grid(starsq,v_xc1Im,vTotIm%pw)
      vTot%pw = vTot%pw + ImagUnit * vTotIm%pw
      call timestop("pw_from_grid")

!      call timestart("finish_pw_grid")
!      CALL finish_pw_grid()
!      call timestop("finish_pw_grid")
   END SUBROUTINE dfpt_vis_xc
END MODULE m_dfpt_vis_xc
