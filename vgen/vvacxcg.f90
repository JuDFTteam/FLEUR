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
  SUBROUTINE vvac_xc(ifftd2,stars,vacuum,noco,oneD,cell,xcpot,input,den, vxc,exc)

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
    USE m_od_mkgxyz3
    USE m_od_mkgz
    USE m_fft2d
    use m_vac_tofrom_grid
    IMPLICIT NONE

    CLASS(t_xcpot),INTENT(IN)    :: xcpot
    TYPE(t_oneD),INTENT(IN)      :: oneD
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
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: rho(:,:),v_xc(:,:),v_x(:,:),e_xc(:,:)

    TYPE(t_gradients)::grad

    !     .. unused input (needed for other noco GGA-implementations) ..

    SELECT TYPE(xcpot)
    TYPE IS (t_xcpot_libxc)
       CALL judft_error("libxc GGA functionals not implemented in film setups")
    END SELECT

    ngrid=vacuum%nvac*(vacuum%nmzxy*ifftd2+vacuum%nmz)

    if (xcpot%needs_grad()) CALL xcpot%alloc_gradients(ngrid,input%jspins,grad)
    allocate(rho(ngrid,input%jspins),v_xc(ngrid,input%jspins),v_x(ngrid,input%jspins))

    call vac_to_grid(xcpot%needs_grad(),ifftd2,input,vacuum,noco,cell,den,stars,rho,grad)


    !         calculate the exchange-correlation potential in  real space
    !
    CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)

    IF (xcpot%needs_grad()) THEN
      SELECT TYPE(xcpot)
      TYPE IS (t_xcpot_libxc)
        !CALL libxc_postprocess_gga_vac(xcpot,stars,cell,v_xc,grad)
      END SELECT
    ENDIF

    call vac_from_grid(stars,vacuum,v_xc,ifftd2,vxc%vacz,vxc%vacxy)


    IF (ALLOCATED(exc%vacz)) THEN
      ALLOCATE ( e_xc(ngrid,1) ); e_xc=0.0
      CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad, mt_call=.False.)
      CALL vac_from_grid(stars,vacuum,e_xc,ifftd2,exc%vacz,exc%vacxy)
    ENDIF


  END SUBROUTINE vvac_xc
END MODULE m_vvac_xc
