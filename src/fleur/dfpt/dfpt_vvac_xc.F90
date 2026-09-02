MODULE m_dfpt_vvac_xc
    use m_juDFT
    private

    public dfpt_vvac_xc
    !-----------------------------------------------------------------------
    !     calculates 2-d star function coefficients of exchange-correlation*
    !     potential in the vacuum regions and adds them to the corresponding
    !     coeffs of the coulomb potential            c.l.fu, r.podloucky   *
    !     for the gradient contribution.   t.a. 1996
    !-----------------------------------------------------------------------
  CONTAINS
    SUBROUTINE dfpt_vvac_xc(ifftd2,stars,starsq, vacuum, noco,cell,den,den1,xcpot,input,vxc)

      USE m_types
      USE m_types_xcpot_libxc
      use m_constants
      use m_vac_tofrom_grid
      USE m_dfpt_gga_kernel
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_stars),INTENT(IN)     :: stars,starsq
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_potden),INTENT(IN)    :: den,den1
      TYPE(t_potden),INTENT(INOUT) :: vxc

      !     ..
      !     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ifftd2

      !     ..
      !     .. Local Scalars ..
      INTEGER :: i,ngrid
      INTEGER :: iSpin,jSpin,nfxc,fxcSpin,n_sigma
      !     ..
      !     .. Local Arrays ..
      REAL, ALLOCATABLE :: rho(:,:), rho1re(:,:),rho1im(:,:)
      REAL, ALLOCATABLE :: rhoDummy(:,:),rhoDummy1(:,:),rhoDummy1Im(:,:)
      REAL, ALLOCATABLE :: v_xc1re(:,:),v_xc1im(:,:),f_xc(:,:)
      !Kernels come out of libxc with the spin-like index first; drivsigmaT is the
      !transpose needed to put drivsigma back on the grid.
      REAL, ALLOCATABLE :: drivsigma(:,:),driv2rho2(:,:),driv2rhosigma(:,:),driv2sigma2(:,:)
      REAL, ALLOCATABLE :: drivsigmaT(:,:),drivsigma1(:,:),drivsigma1Im(:,:)
      COMPLEX, ALLOCATABLE :: drivsigmaVac(:,:,:,:),drivsigma1Vac(:,:,:,:),drivsigma1VacIm(:,:,:,:)
      LOGICAL, ALLOCATABLE :: l_cut(:)
      TYPE(t_gradients)::  gradRho, gradRho1, gradRho1Im
      TYPE(t_gradients)::  gradDrivsigma, gradDrivsigma1, gradDrivsigma1Im
      TYPE(t_potden)   :: vxcIm

      nfxc = 2 * input%jspins - 1
      n_sigma = MERGE(1,3,input%jspins==1)

      ngrid=vacuum%nvac*(vacuum%nmzxy*ifftd2+vacuum%nmz)

      IF (xcpot%needs_grad()) THEN
         CALL xcpot%alloc_gradients(ngrid,input%jspins,gradRho)
         CALL xcpot%alloc_gradients(ngrid,input%jspins,gradRho1)
         !The grid is shorter than the array beyond the warping region
         gradRho1%gr = 0.0
         gradRho1%laplace = 0.0
      END IF
      allocate(rho(ngrid,input%jspins))
      allocate(rho1re(ngrid,input%jspins))
      allocate(rho1im(ngrid,input%jspins))

      ALLOCATE(v_xc1re(ngrid,input%jspins))
      ALLOCATE(v_xc1im(ngrid,input%jspins))

      CALL vxcIm%copyPotDen(vxc)
      CALL vxcIm%resetPotDen()

      rho=0.0
      rho1re=0.0
      rho1im=0.0
      !den1 lives on the q-shifted stars, so its in-plane derivatives carry
      !i(q+G) and the response gradient is complex.
      call timestart("vac_to_grid")
      call vac_to_grid(xcpot%needs_grad(),ifftd2,input%jspins,vacuum,noco%l_noco,cell,den%vac,stars,rho,gradRho)
      IF (xcpot%needs_grad()) THEN
         call vac_to_grid(.TRUE.,ifftd2,input%jspins,vacuum,noco%l_noco,cell,den1%vac,starsq,rho1re,gradRho1,rho1im,gradim=gradRho1Im)
      ELSE
         call vac_to_grid(.FALSE.,ifftd2,input%jspins,vacuum,noco%l_noco,cell,den1%vac,starsq,rho1re,gradRho1,rho1im)
      END IF
      call timestop("vac_to_grid")

      IF (xcpot%needs_grad()) THEN
         !Same cutoff as vvac_xc, so the kernels stay the derivative of the
         !potential it builds; the response is dropped where the ground state is.
         call timestart("apply_vac_cutoffs")
         SELECT TYPE(xcpot)
         TYPE IS (t_xcpot_libxc)
            CALL xcpot%apply_vac_cutoffs(rho,gradRho,l_cut)
         END SELECT
         DO i = 1, ngrid
            IF (l_cut(i)) THEN
               gradRho1%gr(:,i,:) = 0.0
               gradRho1%laplace(i,:) = 0.0
               gradRho1Im%gr(:,i,:) = 0.0
               gradRho1Im%laplace(i,:) = 0.0
            END IF
         END DO
         call timestop("apply_vac_cutoffs")

         ALLOCATE(drivsigma(n_sigma,ngrid),driv2rho2(nfxc,ngrid))
         ALLOCATE(driv2rhosigma(input%jspins*n_sigma,ngrid),driv2sigma2(MERGE(1,6,input%jspins==1),ngrid))
         SELECT TYPE(xcpot)
         TYPE IS (t_xcpot_libxc)
            CALL xcpot%get_fxc_gga(input%jspins,rho,gradRho%sigma,drivsigma,driv2rho2,driv2rhosigma,driv2sigma2)
         END SELECT

         !The cut points carry no density, so libxc returns zero kernels there;
         !this only keeps that explicit for the star expansions below.
         DO i = 1, ngrid
            IF (l_cut(i)) drivsigma(:,i) = 0.0
         END DO

         ALLOCATE(drivsigma1(ngrid,n_sigma),drivsigma1Im(ngrid,n_sigma))
         call timestart("dfpt_gga_assemble")
         CALL dfpt_gga_assemble(drivsigma,driv2rho2,driv2rhosigma,driv2sigma2,gradRho,gradRho1,rho1re,v_xc1re,drivsigma1)
         CALL dfpt_gga_assemble(drivsigma,driv2rho2,driv2rhosigma,driv2sigma2,gradRho,gradRho1Im,rho1im,v_xc1im,drivsigma1Im)
         call timestop("dfpt_gga_assemble")

         !The floored density is a constant, so its response vanishes there.
         DO i = 1, ngrid
            IF (l_cut(i)) THEN
               drivsigma1(i,:) = 0.0
               drivsigma1Im(i,:) = 0.0
               v_xc1re(i,:) = 0.0
               v_xc1im(i,:) = 0.0
            END IF
         END DO

         call timestart("gradient of drivsigma")
         !Ground state drivsigma; unshifted stars and real.
         ALLOCATE(drivsigmaT(ngrid,n_sigma),rhoDummy(ngrid,n_sigma))
         ALLOCATE(drivsigmaVac(vacuum%nmzd,stars%ng2,vacuum%nvac,n_sigma))
         drivsigmaT = TRANSPOSE(drivsigma)
         drivsigmaVac = CMPLX(0.0,0.0)
         CALL vac_from_grid(stars,vacuum,drivsigmaT,ifftd2,drivsigmaVac)
         ALLOCATE(gradDrivsigma%gr(3,ngrid,n_sigma),gradDrivsigma%sigma(n_sigma,ngrid))
         CALL vac_to_grid(.TRUE.,ifftd2,n_sigma,vacuum,.FALSE.,cell,drivsigmaVac,stars,rhoDummy,gradDrivsigma)

         !Response of drivsigma; q-shifted stars, so both channels are needed.
         ALLOCATE(rhoDummy1(ngrid,n_sigma))
         ALLOCATE(drivsigma1Vac(vacuum%nmzd,starsq%ng2,vacuum%nvac,n_sigma))
         ALLOCATE(drivsigma1VacIm(vacuum%nmzd,starsq%ng2,vacuum%nvac,n_sigma))
         drivsigma1Vac = CMPLX(0.0,0.0)
         drivsigma1VacIm = CMPLX(0.0,0.0)
         CALL vac_from_grid(starsq,vacuum,drivsigma1,ifftd2,drivsigma1Vac)
         CALL vac_from_grid(starsq,vacuum,drivsigma1Im,ifftd2,drivsigma1VacIm)
         drivsigma1Vac = drivsigma1Vac + ImagUnit * drivsigma1VacIm
         ALLOCATE(gradDrivsigma1%gr(3,ngrid,n_sigma),gradDrivsigma1%sigma(n_sigma,ngrid))
         CALL vac_to_grid(.TRUE.,ifftd2,n_sigma,vacuum,.FALSE.,cell,drivsigma1Vac,starsq,rhoDummy1,gradDrivsigma1,rhoDummy1Im,gradim=gradDrivsigma1Im)
         call timestop("gradient of drivsigma")

         call timestart("dfpt_gga_grdotgr")
         CALL dfpt_gga_grdotgr(gradRho,gradRho1,gradDrivsigma,gradDrivsigma1,v_xc1re)
         CALL dfpt_gga_grdotgr(gradRho,gradRho1Im,gradDrivsigma,gradDrivsigma1Im,v_xc1im)
         call timestop("dfpt_gga_grdotgr")
      ELSE
         ALLOCATE(f_xc(ngrid,nfxc))
#ifdef CPP_LIBXC
         CALL xcpot%get_fxc_lda(input%jspins, rho, f_xc)
#endif

         v_xc1re = 0.0
         v_xc1im = 0.0
         DO iSpin = 1, input%jspins
             DO jSpin = 1, input%jspins
                 fxcSpin = iSpin + jSpin - 1
                 v_xc1re(:, iSpin) = v_xc1re(:, iSpin) + f_xc(:, fxcSpin) * rho1re(:, jSpin)
                 v_xc1im(:, iSpin) = v_xc1im(:, iSpin) + f_xc(:, fxcSpin) * rho1im(:, jSpin)
             END DO
         END DO
      END IF

      call timestart("vac_from_grid")
      call vac_from_grid(starsq,vacuum,v_xc1re,ifftd2,vxc%vac)
      call vac_from_grid(starsq,vacuum,v_xc1im,ifftd2,vxcIm%vac)
      vxc%vac=vxc%vac + ImagUnit * vxcIm%vac
      call timestop("vac_from_grid")

    END SUBROUTINE dfpt_vvac_xc
  END MODULE m_dfpt_vvac_xc
