MODULE m_vvac_xc
  use m_juDFT
  private
  !These used to be inputs for testing...
  INTEGER,PARAMETER:: fixed_ndvgrd=6
  REAL,PARAMETER   :: fixed_chng=-0.1e-11

  public vvac_xc
  !-----------------------------------------------------------------------
  !     calculates 2-d star function coefficients of exchange-correlation*
  !     potential in the vacuum regions and adds them to the corresponding
  !     coeffs of the coulomb potential            c.l.fu, r.podloucky   *
  !     for the gradient contribution.   t.a. 1996
  !-----------------------------------------------------------------------
CONTAINS
  SUBROUTINE vvac_xc(ifftd2,stars,vacuum,noco ,cell,xcpot,input,den, vxc,exc)

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
     
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_potden),INTENT(IN)    :: den
    TYPE(t_potden),INTENT(INOUT) :: vxc
    TYPE(t_potden),INTENT(INOUT) :: exc
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ifftd2

    !     ..
    !     .. Local Scalars ..
    INTEGER :: js,nt,i,iq,irec2,nmz0,nmzdiff,ivac,ip,ngrid
    REAL    :: rhti,zro,fgz,rhmnv,d_15,bmat1(3,3),rd
    LOGICAL :: l_libxc
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: rho(:,:),v_xc(:,:),v_x(:,:),e_xc(:,:)

    TYPE(t_gradients)::grad

    !     .. unused input (needed for other noco GGA-implementations) ..
    
    l_libxc=.FALSE.

    !SELECT TYPE(xcpot)
    !TYPE IS (t_xcpot_libxc)
    !   IF (xcpot%needs_grad()) THEN
    !      CALL judft_error("libxc GGA functionals not implemented in film setups")
    !   END IF
    !END SELECT

    ngrid=vacuum%nvac*(vacuum%nmzxy*ifftd2+vacuum%nmz)

    if (xcpot%needs_grad()) CALL xcpot%alloc_gradients(ngrid,input%jspins,grad)
    allocate(rho(ngrid,input%jspins),v_xc(ngrid,input%jspins),v_x(ngrid,input%jspins))
    rho=0.0
    call vac_to_grid(xcpot%needs_grad(),ifftd2,input%jspins,vacuum,noco%l_noco,cell,den%vac,stars,rho,grad)

    !         calculate the exchange-correlation potential in  real space
    !
    CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)

    SELECT TYPE(xcpot)
    TYPE IS (t_xcpot_libxc)
       l_libxc=.TRUE.
       IF (xcpot%needs_grad()) THEN
          CALL libxc_postprocess_gga_vac(xcpot,input,cell,stars,vacuum ,v_xc,grad)
          CALL libxc_postprocess_gga_vac(xcpot,input,cell,stars,vacuum ,v_x,grad)
       END IF
    END SELECT

    call vac_from_grid(stars,vacuum,v_xc,ifftd2,vxc%vac)

    !IF (l_libxc.AND.xcpot%needs_grad()) THEN
    !   CALL save_npy('vxc_gga_vac_libxc.npy',v_xc)
    !ELSE IF (l_libxc.AND.(.NOT.xcpot%needs_grad())) THEN
    !  CALL save_npy('vxc_lda_vac_libxc.npy',v_xc)
    !ELSE IF ((.NOT.l_libxc).AND.xcpot%needs_grad()) THEN
    !   CALL save_npy('vxc_gga_vac_inbuild.npy',v_xc)
    !ELSE
    !  CALL save_npy('vxc_lda_vac_inbuild.npy',v_xc)
    !END IF

    IF (ALLOCATED(exc%vac)) THEN
      ALLOCATE ( e_xc(ngrid,1) ); e_xc=0.0
      CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad, mt_call=.False.)
      CALL vac_from_grid(stars,vacuum,e_xc,ifftd2,exc%vac)
    ENDIF


  END SUBROUTINE vvac_xc
END MODULE m_vvac_xc
