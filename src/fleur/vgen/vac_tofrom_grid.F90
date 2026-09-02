!--------------------------------------------------------------------------------
! Copyright (C) 2020 Peter GrüNberg Institut, Forschungszentrum JüLich, Germany
! This File Is Part Of Fleur And Available As Free Software Under The Conditions
! Of The Mit License As Expressed In The License File In More Detail.
!--------------------------------------------------------------------------------
MODULE m_vac_tofrom_grid
      INTEGER,PARAMETER :: fixed_ndvgrd=6

CONTAINS
  subroutine vac_to_grid(dograds,ifftd2,jspins,vacuum,l_noco,cell,vacnew,stars,rho,grad,rhoim,gradim)


    !-----------------------------------------------------------------------
    !     instead of vvacxcor.f: the different exchange-correlation
    !     potentials defined through the key icorr are called through
    !     the driver subroutine vxcallg.f, subroutines vectorized
    !     in case of total = .true. calculates the ex-corr. energy
    !     density through the driver subroutine excallg.f
    !     ** r.pentcheva 08.05.96
    !-----------------------------------------------------------------------

    USE m_juDFT
    USE m_types
    use m_constants
    USE m_grdrsvac
    USE m_grdchlh
    USE m_mkgz
    USE m_mkgxyz3
    !
    !
    USE m_fft2d

    IMPLICIT NONE
    logical,intent(in)           :: dograds
    INTEGER,INTENT(IN)           :: jspins
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    LOGICAL,INTENT(IN)           :: l_noco
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    COMPLEX,INTENT(IN)    :: vacnew(:,:,:,:)
    TYPE(t_gradients),INTENT(INOUT)::grad
    real,intent(OUT)             :: rho(:,:)
    real, optional, allocatable, intent(out) :: rhoim(:,:)
    !Imaginary part of the derivatives; only for a complex (DFPT) field
    TYPE(t_gradients), OPTIONAL, INTENT(INOUT) :: gradim
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ifftd2

    !     ..
    !     .. Local Scalars ..
    INTEGER :: js,nt,i,iq,irec2,nmz0,nmzdiff,ivac,ip,idx,idx1,idx_loc,ipt
    REAL    :: rhti,zro,fgz,rhmnv,d_15,rd,qsq
    LOGICAL :: l_gradim
    !     ..
    !     .. Local Arrays ..
    REAL    :: qvec(3)
    REAL, ALLOCATABLE :: bf2(:)
    REAL, ALLOCATABLE :: rhdx(:,:),rhdy(:,:),rhdz(:,:)
    REAL, ALLOCATABLE :: rhdxx(:,:),rhdyy(:,:),rhtdz(:,:),rhtdzz(:,:)
    REAL, ALLOCATABLE :: rhdzz(:,:),rhdyz(:,:),rhdzx(:,:),rhdxy(:,:)
    REAL, ALLOCATABLE :: rhdxim(:,:),rhdyim(:,:),rhdzim(:,:)
    REAL, ALLOCATABLE :: rhdxxim(:,:),rhdyyim(:,:),rhdzzim(:,:)
    REAL, ALLOCATABLE :: rhdyzim(:,:),rhdzxim(:,:),rhdxyim(:,:)
    REAL, ALLOCATABLE :: rhtdzim(:,:),rhtdzzim(:,:)
    REAL, ALLOCATABLE :: rxydzr(:),rxydzi(:)
    REAL, ALLOCATABLE :: rxydzzr(:),rxydzzi(:),rhtxyr(:),rhtxyi(:)
    REAL, ALLOCATABLE :: rhtxc(:,:)
    COMPLEX, ALLOCATABLE :: rxydz(:,:,:),rxydzz(:,:,:),cqpw(:)
    COMPLEX, ALLOCATABLE :: rdz(:,:,:),rdzz(:,:,:)

    !     ..
    !     for the noco-case only
    REAL :: chdens
    REAL, ALLOCATABLE :: magmom(:,:), dxmagmom(:),ddxmagmom(:,:)
    REAL, ALLOCATABLE :: dymagmom(:),ddymagmom(:,:), dzmagmom(:,:),ddzmagmom(:,:)
    REAL, ALLOCATABLE :: mx(:),my(:)
    !     .. unused input (needed for other noco GGA-implementations) ..



    d_15     = 1.e-15
    zro      = 0.0
    nt       = ifftd2
    idx=1

    rho = 0.0

    !stars%center is zero except for the q-shifted stars used in DFPT; q is
    !in-plane there, so only the x and y derivatives pick up a shift.
    qvec = matmul(stars%center,cell%bmat)
    qsq = dot_product(qvec(:2),qvec(:2))

    l_gradim = PRESENT(gradim)
    IF (l_gradim.AND.l_noco) CALL judft_error("vac_to_grid: complex gradients are not available for noco",calledby="vac_tofrom_grid.F90")
    IF (l_gradim.AND..NOT.PRESENT(rhoim)) CALL judft_error("vac_to_grid: gradim needs rhoim",calledby="vac_tofrom_grid.F90")

    ALLOCATE ( bf2(ifftd2) )
    IF (PRESENT(rhoim)) THEN
      ALLOCATE(rhoim,mold=rho)
      rhoim=0.0
    ENDIF
    IF (l_gradim) THEN
      IF (ALLOCATED(gradim%gr)) DEALLOCATE(gradim%gr)
      IF (ALLOCATED(gradim%laplace)) DEALLOCATE(gradim%laplace)
      ALLOCATE(gradim%gr(3,SIZE(rho,1),jspins),gradim%laplace(SIZE(rho,1),jspins))
      gradim%gr = 0.0
      gradim%laplace = 0.0
    ENDIF
    !The grid is longer than the part the vacua actually fill, so the tail would
    !otherwise reach the functionals uninitialised.
    IF (ALLOCATED(grad%gr)) grad%gr = 0.0
    IF (ALLOCATED(grad%laplace)) grad%laplace = 0.0
    IF (ALLOCATED(grad%sigma)) grad%sigma = 0.0
    WRITE (oUnit,'(/'' ifftd2,vacuum%nmz='',2i7)') ifftd2,vacuum%nmz
    WRITE (oUnit,'('' 9990nmzxy='',2i5)') vacuum%nmzxy

    ALLOCATE ( rxydz(vacuum%nmzxy,stars%ng2,jspins),rxydzz(vacuum%nmzxyd,stars%ng2,jspins) )
    ALLOCATE ( rhtdz(vacuum%nmzd,jspins),rhtdzz(vacuum%nmzd,jspins) )
    ALLOCATE ( rdz(vacuum%nmzd,stars%ng2,jspins),rdzz(vacuum%nmzd,stars%ng2,jspins))
    IF (l_gradim) ALLOCATE ( rhtdzim(vacuum%nmzd,jspins),rhtdzzim(vacuum%nmzd,jspins) )
    !ALLOCATE ( fgxy(stars%ng2-1) )
	 rxydz = CMPLX(0.0,0.0)
	 rxydzz= CMPLX(0.0,0.0)

    IF (l_noco) THEN
      ALLOCATE ( magmom(0:ifftd2-1,vacuum%nmzxy) )
      ALLOCATE ( dzmagmom(0:ifftd2-1,vacuum%nmzxy) )
      ALLOCATE ( ddzmagmom(0:ifftd2-1,vacuum%nmzxy) )
      ALLOCATE ( mx(0:ifftd2-1),my(0:ifftd2-1) )
    ENDIF
    IF ( l_noco .OR. dograds ) THEN
      ALLOCATE ( rhtxyr(vacuum%nmzxy)  )
      ALLOCATE ( rxydzr(vacuum%nmzxy),rxydzzr(vacuum%nmzxy) )
    ENDIF
    IF (dograds) THEN
      ALLOCATE ( rhtxyi(vacuum%nmzxy) )
      ALLOCATE ( rxydzi(vacuum%nmzxy) )
      ALLOCATE ( rxydzzi(vacuum%nmzxy) )
    ENDIF
    DO ivac=1,vacuum%nvac

       ! the charge density in vacuum is expanded in 2-dim stars on a mesh
       ! in z-direction. the g||.ne.zero-components expand from 1 to nmzxy
       ! the g||.eq.zero-components expand from 1 to nmz
       ! first we calculate vxc in the warping region


          ! Transform charge and magnetization to real-space.
          ! In the collinear case that is done later within
          ! another loop over the vacuum-layers in order to
          ! save memory.

          idx1=idx
          !idx1=(ivac-1)* ( vacuum%nmzxy * ifftd2 + nmzdiff ) + 1
          DO ip=1,vacuum%nmzxy
            DO js=1,jspins
              IF (.NOT.PRESENT(rhoim)) THEN
                 CALL fft2d(stars, rho(idx1:idx1+ifftd2-1,js),bf2, vacnew(ip,:,ivac,js),+1)
              ELSE
                 CALL fft2d(stars, rho(idx1:idx1+ifftd2-1,js),rhoim(idx1:idx1+ifftd2-1,js), vacnew(ip,:,ivac,js),+1)
              END IF
            END DO
            IF (l_noco) THEN
              CALL fft2d(stars, mx,my, vacnew(ip,:,ivac,3),+1)

              DO i=0,ifftd2-1
                magmom(i,ip)= mx(i)**2 + my(i)**2 + ((rho(i+idx1,1)-rho(i+idx1,2))/2.)**2
                magmom(i,ip)= SQRT(magmom(i,ip))
                chdens= rho(i+idx1,1)/2.+rho(i+idx1,2)/2.
                rho(i+idx1,1)= chdens + magmom(i,ip)
                rho(i+idx1,2)= chdens - magmom(i,ip)
              END DO
            ENDIF
            idx1=idx1+ifftd2
          END DO ! ip=1,vacuum%nmzxy

       !      ENDDO    ! ivac
       !      DO ivac = 1,nvac

       IF (dograds) THEN
          DO js=1,jspins
             !
             ! calculate first (rhtdz) & second (rhtdzz) derivative of vacz(1:nmz)
             CALL grdchlh(vacuum%delz,REAL(vacnew(1:vacuum%nmz,1,ivac,js)),&
                  rhtdz(1:,js),rhtdzz(1:,js))
				 rdz(:,1,js) = rhtdz(1:,js)
				 rdzz(:,1,js) = rhtdzz(1:,js)
             IF (l_gradim) THEN
                CALL grdchlh(vacuum%delz,AIMAG(vacnew(1:vacuum%nmz,1,ivac,js)),&
                     rhtdzim(1:,js),rhtdzzim(1:,js))
                rdz(:,1,js) = CMPLX(rhtdz(1:,js),rhtdzim(1:,js))
                rdzz(:,1,js) = CMPLX(rhtdzz(1:,js),rhtdzzim(1:,js))
             END IF
             DO iq = 1, stars%ng2-1
                !
                ! calculate first (rxydz) & second (rxydzz) derivative of vacxy:
                !
                DO ip=1,vacuum%nmzxy
                   rhtxyr(ip)=real(vacnew(ip,iq+1,ivac,js))
                ENDDO
                CALL grdchlh(vacuum%delz,rhtxyr(:vacuum%nmzxy), rxydzr,rxydzzr)

                DO ip=1,vacuum%nmzxy
                   rhtxyi(ip)=aimag(vacnew(ip,iq+1,ivac,js))
                ENDDO
                CALL grdchlh(vacuum%delz,rhtxyi(:vacuum%nmzxy), rxydzi,rxydzzi)

                DO ip=1,vacuum%nmzxy
                   rxydz(ip,iq+1,js)=cmplx(rxydzr(ip),rxydzi(ip))
                   rxydzz(ip,iq+1,js)=cmplx(rxydzzr(ip),rxydzzi(ip))
                ENDDO

             ENDDO ! loop over 2D stars (iq)
				 rdz(:vacuum%nmzxy,2:,js) = rxydz(:vacuum%nmzxy,2:,js)
				 rdzz(:vacuum%nmzxy,2:,js) = rxydzz(:vacuum%nmzxy,2:,js)

          ENDDO ! jspins

          IF (l_noco) THEN
             !  calculate  dzmagmom = d magmom / d z  and ddzmagmom= d dmagmom / d z

             DO i=0,ifftd2-1
                DO ip=1,vacuum%nmzxy
                   rhtxyr(ip)=magmom(i,ip)
                ENDDO
                CALL grdchlh(vacuum%delz,rhtxyr(1:vacuum%nmzxy), rxydzr,rxydzzr)
                DO ip=1,vacuum%nmzxy
                   dzmagmom(i,ip)= rxydzr(ip)
                   ddzmagmom(i,ip)= rxydzzr(ip)
                ENDDO
             END DO
          END IF ! l_noco

       ENDIF   ! xcpot%igrd.GT.0

       !       WRITE(oUnit,'('' 9990nmzxy='',2i5)') nmzxy

       CALL timestart("warp")
       rd = 0.0
       !$OMP PARALLEL DEFAULT(none) &
       !$OMP SHARED(vacuum,dograds,jspins,stars,ivac,zro,cell,magmom,vacnew) &
       !$OMP SHARED(rhtdz,rhtdzz,rdz,rdzz,rxydz,rxydzz,l_noco,dzmagmom,ddzmagmom,idx) &
       !$OMP SHARED(ifftd2,rho,rhoim,grad,gradim,l_gradim) &
       !$OMP PRIVATE(ip,js,iq,cqpw,rhti,rhdx,rhdy,rhdz,rhdxx,rhdyy,rhdzz) &
       !$OMP PRIVATE(rhdxy,rhdzx,rhdyz,dxmagmom,dymagmom,ddxmagmom,ddymagmom) &
       !$OMP PRIVATE(rhdxim,rhdyim,rhdzim,rhdxxim,rhdyyim,rhdzzim) &
       !$OMP PRIVATE(rhdxyim,rhdzxim,rhdyzim) &
       !$OMP PRIVATE(chdens,idx_loc)
       ALLOCATE ( rhdx(0:ifftd2-1,jspins),rhdy(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdz(0:ifftd2-1,jspins),rhdxx(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdyy(0:ifftd2-1,jspins),rhdzz(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdyz(0:ifftd2-1,jspins),rhdzx(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdxy(0:ifftd2-1,jspins))
       ALLOCATE ( rhdxim(0:ifftd2-1,jspins),rhdyim(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdzim(0:ifftd2-1,jspins),rhdxxim(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdyyim(0:ifftd2-1,jspins),rhdzzim(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdyzim(0:ifftd2-1,jspins),rhdzxim(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdxyim(0:ifftd2-1,jspins))
       ALLOCATE ( cqpw(stars%ng2))
       IF (l_noco) THEN
          ALLOCATE ( dxmagmom(0:ifftd2-1),dymagmom(0:ifftd2-1) )
          ALLOCATE ( ddxmagmom(0:ifftd2-1,2),ddymagmom(0:ifftd2-1,2) )
       ENDIF
       !$OMP DO
       DO ip = 1,vacuum%nmzxy
          ! loop over warping region


          IF (dograds) THEN
             ! calculate derivatives with respect to x,y in g-space
             ! and transform them to real-space.

             DO js = 1,jspins

                !The g||=0 star is kept: with the DFPT q-shift its in-plane
                !derivative is i*q, and for stars%center=0 the factor vanishes anyway.
                cqpw(1:stars%ng2)=vacnew(ip,1:stars%ng2,ivac,js)

                CALL fft2d(stars, rhdx(0,js),rhdxim(0,js), cqpw, +1,firstderiv=[1.,0.,0.],cell=cell)
                CALL fft2d(stars, rhdy(0,js),rhdyim(0,js), cqpw, +1,firstderiv=[0.,1.,0.],cell=cell) ! dn/dy =  FFT(0,i*gy*vacxy)&
                CALL fft2d(stars, rhdz(0,js),rhdzim(0,js), rdz(ip,:,js), +1) ! dn/dz = FFT(rhtdz,rxydz)&

                CALL fft2d(stars, rhdxx(0,js),rhdxxim(0,js), cqpw, +1,firstderiv=[1.0,0.,0.],secondderiv=[1.0,0.,0.],cell=cell) ! d2n/dx2 = FFT(0,-gx^2*vacxy)&
                CALL fft2d(stars, rhdyy(0,js),rhdyyim(0,js), cqpw, +1,firstderiv=[0.,1.0,0.],secondderiv=[0.,1.0,0.],cell=cell) ! d2n/dy2 = FFT(0,-gy^2*vacxy)&
                CALL fft2d(stars, rhdzz(0,js),rhdzzim(0,js), rdzz(ip,:,js), +1) ! d2n/dz2 = FFT(rhtdzz,rxydzz)&
                CALL fft2d(stars, rhdxy(0,js),rhdxyim(0,js), cqpw, +1,firstderiv=[0.,1.0,0.],secondderiv=[1.,0.0,0.],cell=cell) ! d2n/dxy = FFT(0,-gx*gy*vacxy)&


                cqpw(1:stars%ng2)=rdz(ip,1:stars%ng2,js)

                CALL fft2d(stars, rhdyz(0,js),rhdyzim(0,js), cqpw, +1,firstderiv=[0.,1.0,0.],cell=cell) ! d2n/dyz = FFT(0,i*gy*rxydz)&
                CALL fft2d(stars, rhdzx(0,js),rhdzxim(0,js), cqpw, +1,firstderiv=[1.,0.0,0.],cell=cell) ! d2n/dzx = FFT(0,i*gx*rxydz)&


             END DO ! js=1,jspins


             IF (l_noco) THEN
                ! ! In non-collinear calculations the derivatives of |m| are calculated
                ! ! in real-space. The derivatives of the charge density, that are
                ! ! already calculated in g-space, will be used.

                CALL grdrsvac(magmom(0,ip),cell%bmat,3*stars%mx1,3*stars%mx2,fixed_ndvgrd, dxmagmom,dymagmom)
                DO i=0,ifftd2-1
                   chdens= rhdx(i,1)/2.+rhdx(i,2)/2.
                   rhdx(i,1)= chdens + dxmagmom(i)
                   rhdx(i,2)= chdens - dxmagmom(i)
                   chdens= rhdy(i,1)/2.+rhdy(i,2)/2.
                   rhdy(i,1)= chdens + dymagmom(i)
                   rhdy(i,2)= chdens - dymagmom(i)
                   chdens= rhdz(i,1)/2.+rhdz(i,2)/2.
                   rhdz(i,1)= chdens + dzmagmom(i,ip)
                   rhdz(i,2)= chdens - dzmagmom(i,ip)
                END DO

                CALL grdrsvac(dxmagmom,cell%bmat,3*stars%mx1,3*stars%mx2,fixed_ndvgrd, &
                     ddxmagmom(0,1),ddymagmom(0,1))
                CALL grdrsvac(&
                     dymagmom,cell%bmat,3*stars%mx1,3*stars%mx2,fixed_ndvgrd,ddxmagmom(0,2),ddymagmom(0,2))
                DO i=0,ifftd2-1
                   chdens= rhdxx(i,1)/2.+rhdxx(i,2)/2.
                   rhdxx(i,1)= chdens + ddxmagmom(i,1)
                   rhdxx(i,2)= chdens - ddxmagmom(i,1)
                   chdens= rhdyy(i,1)/2.+rhdyy(i,2)/2.
                   rhdyy(i,1)= chdens + ddymagmom(i,2)
                   rhdyy(i,2)= chdens - ddymagmom(i,2)
                   chdens= rhdxy(i,1)/2.+rhdxy(i,2)/2.
                   rhdxy(i,1)= chdens + ( ddxmagmom(i,2) + ddymagmom(i,1) )/2.
                   rhdxy(i,2)= chdens - ( ddxmagmom(i,2) + ddymagmom(i,1) )/2.
                END DO
                CALL grdrsvac(dzmagmom(0,ip),cell%bmat,3*stars%mx1,3*stars%mx2,fixed_ndvgrd, &
                     ddxmagmom(0,1),ddymagmom(0,1))
                DO i=0,ifftd2-1
                   chdens= rhdzx(i,1)/2.+rhdzx(i,2)/2.
                   rhdzx(i,1)= chdens + ddxmagmom(i,1)
                   rhdzx(i,2)= chdens - ddxmagmom(i,1)
                   chdens= rhdyz(i,1)/2.+rhdyz(i,2)/2.
                   rhdyz(i,1)= chdens + ddymagmom(i,1)
                   rhdyz(i,2)= chdens - ddymagmom(i,1)
                   chdens= rhdzz(i,1)/2.+rhdzz(i,2)/2.
                   rhdzz(i,1)= chdens + ddzmagmom(i,ip)
                   rhdzz(i,2)= chdens - ddzmagmom(i,ip)
                END DO

             END IF ! l_noco
!        
             idx_loc = idx + (ip-1)* ifftd2
             CALL mkgxyz3(rho(idx_loc:idx_loc+ifftd2-1,:),rhdx,rhdy, rhdz,rhdxx,rhdyy,rhdzz,rhdyz,rhdzx,rhdxy, idx_loc-1,grad)
             IF (l_gradim) CALL mkgxyz3(0*rhdxim,rhdxim,rhdyim,rhdzim,rhdxxim,rhdyyim,rhdzzim,rhdyzim,rhdzxim,rhdxyim, idx_loc-1,gradim)
!
          END IF ! vxc_is_gga
          !
          ! set minimal value of af2 to 1.0e-13
          !

!          rho=max(rho,1e-13)

       END DO ! ip=1,vacuum%nmzxy
       !$OMP END DO
       DEALLOCATE ( rhdx,rhdy )
       DEALLOCATE ( rhdz,rhdxx )
       DEALLOCATE ( rhdyy,rhdzz )
       DEALLOCATE ( rhdyz,rhdzx )
       DEALLOCATE ( rhdxy,cqpw )
       DEALLOCATE ( rhdxim,rhdyim,rhdzim,rhdxxim )
       DEALLOCATE ( rhdyyim,rhdzzim,rhdyzim,rhdzxim,rhdxyim )
       IF (l_noco) THEN
          DEALLOCATE ( dxmagmom,dymagmom )
          DEALLOCATE ( ddxmagmom,ddymagmom )
       ENDIF
       !$OMP END PARALLEL 
       idx = idx + vacuum%nmzxy * ifftd2
       CALL timestop("warp")

       ! now treat the non-warping region


       ! The non-warping region runs from nmzxy+1 to nmz.
       ! The values from nmz0 to nmzxy are taken into account in order
       ! to get the real-space derivative smooth around nmzxy+1.
       nmz0= max(1,vacuum%nmzxy+1+(fixed_ndvgrd/2)-fixed_ndvgrd)
       nmzdiff = vacuum%nmz - nmz0+1
       !       WRITE(oUnit,'(/'' 9992excz''/(8f15.7))') (excz(ip,1),ip=1,nmz)
       WRITE(oUnit,'(/'' 9992nmzdiff='',i5)') nmzdiff


       !idx = (ivac-1)* ( vacuum%nmzxy * ifftd2 + nmzdiff ) + ip*ifftd2 + 1
       DO ip=nmz0,vacuum%nmz
          IF (.not. l_noco) THEN
             DO js=1,jspins
                !rho(idx+ip-nmz0,js)= REAL(vacz(ip,ivac,js))
                rho(idx+ip-nmz0,js)= REAL(vacnew(ip,1,ivac,js))
                IF (PRESENT(rhoim)) rhoim(idx+ip-nmz0,js)= AIMAG(vacnew(ip,1,ivac,js))
             END DO
          ELSE
             !mx(0) = REAL(vacz(ip,ivac,3))
             mx(0) = REAL(vacnew(ip,1,ivac,3))
             !my(0) = AIMAG(vacz(ip,ivac,3))
             my(0) = AIMAG(vacnew(ip,1,ivac,3))
             !chdens= (REAL(vacz(ip,ivac,1))+REAL(vacz(ip,ivac,2)))/2
             chdens= (REAL(vacnew(ip,1,ivac,1))+REAL(vacnew(ip,1,ivac,2)))/2
             !magmom(0,1)= mx(0)**2 + my(0)**2 + ((REAL(vacz(ip,ivac,1))-REAL(vacz(ip,ivac,2)))/2.)**2
             magmom(0,1)= mx(0)**2 + my(0)**2 + ((REAL(vacnew(ip,1,ivac,1))-REAL(vacnew(ip,1,ivac,2)))/2.)**2
             magmom(0,1)= SQRT(magmom(0,1))
             rho(idx+ip-nmz0,1)= chdens + magmom(0,1)
             rho(idx+ip-nmz0,2)= chdens - magmom(0,1)
          END IF
       END DO
       IF (dograds)  THEN
         IF (l_noco) THEN
           DO js=1,jspins
             CALL grdchlh(vacuum%delz,rho(idx:idx+nmzdiff-1,js),rhtdz(nmz0:,js),rhtdzz(nmz0:,js))
           END DO

         END IF

       !       calculate the quantities such as abs(grad(rho)),.. used in
       !c      evaluating the gradient contributions to potential and
       !c      energy.



             CALL mkgz(nmzdiff,jspins, rho(nmz0:,1),rho(nmz0:,jspins),&
             rhtdz(nmz0:,1),rhtdz(nmz0:,jspins),rhtdzz(nmz0:,1),&
                  rhtdzz(nmz0:,jspins),idx-1,grad)

             !mkgz only fills the contracted gradients, so gr and laplace are done
             !here. Only the g||=0 coefficient survives beyond nmzxy, hence the
             !in-plane derivatives reduce to the DFPT q-shift.
             DO js = 1,jspins
                DO i = 1,nmzdiff
                   ipt = idx+i-1
                   IF (ALLOCATED(grad%gr)) THEN
                      grad%gr(:,ipt,js) = [0.0,0.0,rhtdz(nmz0+i-1,js)]
                      IF (l_gradim) grad%gr(:2,ipt,js) = qvec(:2)*rhoim(ipt,js)
                   END IF
                   IF (ALLOCATED(grad%laplace)) grad%laplace(ipt,js) = rhtdzz(nmz0+i-1,js) - qsq*rho(ipt,js)
                   IF (l_gradim) THEN
                      gradim%gr(:,ipt,js) = [-qvec(1)*rho(ipt,js),-qvec(2)*rho(ipt,js),rhtdzim(nmz0+i-1,js)]
                      gradim%laplace(ipt,js) = rhtdzzim(nmz0+i-1,js) - qsq*rhoim(ipt,js)
                   END IF
                END DO
             END DO

       ENDIF

       !       calculate vxc for z now beyond warping region
       idx=idx+nmzdiff
    ENDDO    ! loop over vacua (ivac)



  END SUBROUTINE vac_to_grid

  subroutine vac_from_grid(stars,vacuum,v_xc,ifft2d,vac)
    use m_types_stars
    use m_types_vacuum
    use m_fft2d
    type(t_stars),intent(in)  :: stars
    type(t_vacuum),intent(in) :: vacuum

    real, INTENT(IN)          :: v_xc(:,:)
    INTEGER,INTENT(IN)        :: ifft2d
    complex,intent(INOUT)     :: vac(:,:,:,:)

    COMPLEX, ALLOCATABLE    :: fg(:)
    REAL, ALLOCATABLE       :: bf2(:)
    INTEGER                 :: js,irec2,idx,ivac,ip,nmz0

    ALLOCATE ( fg(stars%ng2),bf2(ifft2d) )

    DO js = 1,size(v_xc,2)
      idx=1
      DO ivac=1,vacuum%nvac
        DO ip=1,vacuum%nmzxy
          !
          !           ----> 2-d back fft to g space
          !
          bf2=0.0
          CALL fft2d(stars, v_xc(idx:idx-1+ifft2d,js),bf2, fg, -1)
          idx=idx+ifft2d
          !            ----> and add vxc to coulomb potential
          !                  the g||.eq.zero component is added to vxc%vacz
          !
          vac(ip,1,ivac,js) = fg(1) + vac(ip,1,ivac,js)
          !
          !            the g||.ne.zero components are added to vxc%vacxy
          !
          DO irec2 = 2,stars%ng2
            vac(ip,irec2,ivac,js)=vac(ip,irec2,ivac,js)+fg(irec2)
          ENDDO
        enddo

        nmz0= max(1,vacuum%nmzxy+1+(fixed_ndvgrd/2)-fixed_ndvgrd)

        DO ip = nmz0,vacuum%nmz
          if (ip>vacuum%nmzxy) vac(ip,1,ivac,js) = vac(ip,1,ivac,js) + v_xc(idx,js)
          idx=idx+1
        ENDDO
      END DO
    ENDDO

  end subroutine
end module
