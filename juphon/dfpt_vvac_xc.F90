MODULE m_dfpt_vvac_xc
    use m_juDFT
    private
    !These used to be inputs for testing...
    INTEGER,PARAMETER:: fixed_ndvgrd=6
    REAL,PARAMETER   :: fixed_chng=-0.1e-11
  
    public dfpt_vvac_xc
    !-----------------------------------------------------------------------
    !     calculates 2-d star function coefficients of exchange-correlation*
    !     potential in the vacuum regions and adds them to the corresponding
    !     coeffs of the coulomb potential            c.l.fu, r.podloucky   *
    !     for the gradient contribution.   t.a. 1996
    !-----------------------------------------------------------------------
  CONTAINS
    SUBROUTINE dfpt_vvac_xc(ifftd2,stars,starsq, vacuum, noco,cell,den,den1,xcpot,input,vxc)
  
      !-----------------------------------------------------------------------
      !     instead of vvacxcor.f: the different exchange-correlation
      !     potentials defined through the key icorr are called through
      !     the driver subroutine vxcallg.f, subroutines vectorized
      !     in case of total = .true. calculates the ex-corr. energy
      !     density through the driver subroutine excallg.f
      !     ** r.pentcheva 08.05.96
      !-----------------------------------------------------------------------
  
      USE m_types
      USE m_types_xcpot_libxc
      use m_constants
      USE m_grdrsvac
      USE m_grdchlh
      USE m_mkgz
      USE m_mkgxyz3
      ! 
      ! 
      USE m_fft2d
      use m_vac_tofrom_grid
      USE m_libxc_postprocess_gga
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
      INTEGER :: js,nt,i,iq,irec2,nmz0,nmzdiff,ivac,ip,ngrid
      INTEGER :: iSpin,jSpin,nfxc,fxcSpin
      REAL    :: rhti,zro,fgz,rhmnv,d_15,bmat1(3,3),rd
      LOGICAL :: l_libxc
      !     ..
      !     .. Local Arrays ..
      REAL, ALLOCATABLE :: rho(:,:),v_xc(:,:),v_x(:,:),e_xc(:,:), rho1re(:,:),rho1im(:,:)
      REAL, ALLOCATABLE :: v_xc1re(:,:),v_xc1im(:,:),f_xc(:,:)
      TYPE(t_gradients)::  grad, grad1 !TODO:     not sure if we need grad1
      TYPE(t_potden)   :: vxcIm
      !     .. unused input (needed for other noco GGA-implementations) ..
      
      l_libxc=.FALSE.
      nfxc = 2 * input%jspins - 1  

      !SELECT TYPE(xcpot)
      !TYPE IS (t_xcpot_libxc)
      !   IF (xcpot%needs_grad()) THEN
      !      CALL judft_error("libxc GGA functionals not implemented in film setups")
      !   END IF
      !END SELECT
  
      ngrid=vacuum%nvac*(vacuum%nmzxy*ifftd2+vacuum%nmz)
  
      ALLOCATE(f_xc(SIZE(rho,1),nfxc))
      ALLOCATE(v_xc1re,mold=rho)
      ALLOCATE(v_xc1im,mold=rho)

      if (xcpot%needs_grad()) CALL xcpot%alloc_gradients(ngrid,input%jspins,grad)
      allocate(rho(ngrid,input%jspins),v_xc(ngrid,input%jspins),v_x(ngrid,input%jspins))
      allocate(rhoim(ngrid,input%jspins))
      rho=0.0
      rho1re=0.0
      rho1im=0.0
      !call vac_to_grid(xcpot%needs_grad(),ifftd2,input%jspins,vacuum,noco%l_noco,cell,den%vacxy(:,:,:,:),den%vacz,stars,rho,grad)
      call timestart("vac_to_grid")
      call vac_to_grid(xcpot%needs_grad(),ifftd2,input%jspins,vacuum,noco%l_noco,cell,den%vac(:vacuum%nmzxyd,2:,:,:),REAL(den%vac(:,1,:,:)),den%vac,stars,rho,grad)
      call vac_to_grid(xcpot%needs_grad(),ifftd2,input%jspins,vacuum,noco%l_noco,cell,den1%vac(:vacuum%nmzxyd,2:,:,:),REAL(den1%vac(:,1,:,:)),den1%vac,starsq,rho1re,grad1,rho1im)
      call timestop("vac_to_grid")
      !         calculate the exchange-correlation potential in  real space

#ifdef CPP_LIBXC 
      CALL xcpot%get_fxc(input%jspins, rho, f_xc)
#endif
    
      SELECT TYPE(xcpot)
      TYPE IS (t_xcpot_libxc)
         l_libxc=.TRUE.
         IF (xcpot%needs_grad()) THEN
            CALL judft_error("GGA not yet implemented",calledby ="dfpt_vvac_xc")
            CALL libxc_postprocess_gga_vac(xcpot,input,cell,stars,vacuum ,v_xc,grad)
            CALL libxc_postprocess_gga_vac(xcpot,input,cell,stars,vacuum ,v_x,grad)
         END IF
      END SELECT

      v_xc1re = 0.0
      v_xc1im = 0.0
      DO iSpin = 1, input%jspins
          DO jSpin = 1, input%jspins
              fxcSpin = iSpin + jSpin - 1
              v_xc1re(:, iSpin) = v_xc1re(:, iSpin) + f_xc(:, fxcSpin) * rho1re(:, jSpin)
              v_xc1im(:, iSpin) = v_xc1im(:, iSpin) + f_xc(:, fxcSpin) * rho1im(:, jSpin)
          END DO
      END DO      
      call timestart("vac_from_grid")
      call vac_from_grid(starsq,vacuum,v_xc1re,ifftd2,vxc%vac)
      call vac_from_grid(starsq,vacuum,v_xc1im,ifftd2,vxcIm%vac)
      vxc%vac=vxc%vac + ImagUnit * vxcIm%vac
      call timestop("vac_from_grid")



    END SUBROUTINE dfpt_vvac_xc
  END MODULE m_dfpt_vvac_xc
  