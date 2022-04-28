!--------------------------------------------------------------------------------
! Copyright (C) 2020 Peter GrüNberg Institut, Forschungszentrum JüLich, Germany
! This File Is Part Of Fleur And Available As Free Software Under The Conditions
! Of The Mit License As Expressed In The License File In More Detail.
!--------------------------------------------------------------------------------
MODULE m_vac_tofrom_grid
      INTEGER,PARAMETER :: fixed_ndvgrd=6

CONTAINS
  subroutine vac_to_grid(dograds,ifftd2,jspins,vacuum,l_noco,cell,vacxy,vacz,stars,rho,grad)


    !-----------------------------------------------------------------------
    !     instead of vvacxcor.f: the different exchange-correlation
    !     potentials defined through the key icorr are called through
    !     the driver subroutine vxcallg.f, subroutines vectorized
    !     in case of total = .true. calculates the ex-corr. energy
    !     density through the driver subroutine excallg.f
    !     ** r.pentcheva 08.05.96
    !-----------------------------------------------------------------------

    USE m_types
    use m_constants
    USE m_grdrsvac
    USE m_grdchlh
    USE m_mkgz
    USE m_mkgxyz3
    !USE m_od_mkgxyz3
    !USE m_od_mkgz
    USE m_fft2d

    IMPLICIT NONE
    logical,intent(in)           :: dograds
    INTEGER,INTENT(IN)           :: jspins
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    LOGICAL,INTENT(IN)           :: l_noco
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    COMPLEX,INTENT(IN)    :: vacxy(:,:,:,:)
    REAL,INTENT(IN)    :: vacz(:,:,:)
    TYPE(t_gradients),INTENT(INOUT)::grad
    real,intent(OUT)             :: rho(:,:)
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ifftd2

    !     ..
    !     .. Local Scalars ..
    INTEGER :: js,nt,i,iq,irec2,nmz0,nmzdiff,ivac,ip,idx,idx1,idx_loc
    REAL    :: rhti,zro,fgz,rhmnv,d_15,rd
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: bf2(:)
    REAL, ALLOCATABLE :: rhdx(:,:),rhdy(:,:),rhdz(:,:)
    REAL, ALLOCATABLE :: rhdxx(:,:),rhdyy(:,:),rhtdz(:,:),rhtdzz(:,:)
    REAL, ALLOCATABLE :: rhdzz(:,:),rhdyz(:,:),rhdzx(:,:),rhdxy(:,:)
    REAL, ALLOCATABLE :: rxydzr(:),rxydzi(:)
    REAL, ALLOCATABLE :: rxydzzr(:),rxydzzi(:),rhtxyr(:),rhtxyi(:)
    REAL, ALLOCATABLE :: rhtxc(:,:)
    COMPLEX, ALLOCATABLE :: rxydz(:,:,:),rxydzz(:,:,:),cqpw(:)

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

    ALLOCATE ( bf2(ifftd2) )

    WRITE (oUnit,'(/'' ifftd2,vacuum%nmz='',2i7)') ifftd2,vacuum%nmz
    WRITE (oUnit,'('' 9990nmzxy='',2i5)') vacuum%nmzxy

    ALLOCATE ( rxydz(vacuum%nmzxy,stars%ng2-1,jspins),rxydzz(vacuum%nmzxyd,stars%ng2-1,jspins) )
    ALLOCATE ( rhtdz(vacuum%nmzd,jspins),rhtdzz(vacuum%nmzd,jspins) )
    !ALLOCATE ( fgxy(stars%ng2-1) )


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
              CALL fft2d(stars, rho(idx1:idx1+9*stars%mx1*stars%mx2-1,js),bf2, vacz(ip,ivac,js),0.,&
              vacxy(ip,:,ivac,js),+1)
            END DO
            IF (l_noco) THEN
              CALL fft2d(stars, mx,my, vacz(ip,ivac,3),vacz(ip,ivac,4), &
              vacxy(ip,:,ivac,3),+1)

              DO i=0,9*stars%mx1*stars%mx2-1
                magmom(i,ip)= mx(i)**2 + my(i)**2 + ((rho(i+idx1,1)-rho(i+idx1,2))/2.)**2
                magmom(i,ip)= SQRT(magmom(i,ip))
                chdens= rho(i+idx1,1)/2.+rho(i+idx1,2)/2.
                rho(i+idx1,1)= chdens + magmom(i,ip)
                rho(i+idx1,2)= chdens - magmom(i,ip)
              END DO
            ENDIF
            idx1=idx1+9*stars%mx1*stars%mx2
          END DO ! ip=1,vacuum%nmzxy

       !      ENDDO    ! ivac
       !      DO ivac = 1,nvac

       IF (dograds) THEN
          DO js=1,jspins
             !
             ! calculate first (rhtdz) & second (rhtdzz) derivative of vacz(1:nmz)
             !
             CALL grdchlh(vacuum%delz,vacz(1:vacuum%nmz,ivac,js),&
                  rhtdz(1:,js),rhtdzz(1:,js))

             DO iq = 1, stars%ng2-1
                !
                ! calculate first (rxydz) & second (rxydzz) derivative of vacxy:
                !
                DO ip=1,vacuum%nmzxy
                   rhtxyr(ip)=vacxy(ip,iq,ivac,js)
                ENDDO
                CALL grdchlh(vacuum%delz,rhtxyr(:vacuum%nmzxy), rxydzr,rxydzzr)

                DO ip=1,vacuum%nmzxy
                   rhtxyi(ip)=aimag(vacxy(ip,iq,ivac,js))
                ENDDO
                CALL grdchlh(vacuum%delz,rhtxyi(:vacuum%nmzxy), rxydzi,rxydzzi)

                DO ip=1,vacuum%nmzxy
                   rxydz(ip,iq,js)=cmplx(rxydzr(ip),rxydzi(ip))
                   rxydzz(ip,iq,js)=cmplx(rxydzzr(ip),rxydzzi(ip))
                ENDDO

             ENDDO ! loop over 2D stars (iq)


          ENDDO ! jspins

          IF (l_noco) THEN
             !  calculate  dzmagmom = d magmom / d z  and ddzmagmom= d dmagmom / d z

             DO i=0,9*stars%mx1*stars%mx2-1
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
       !$OMP SHARED(vacuum,dograds,jspins,stars,ivac,zro,cell,magmom,vacxy) &
       !$OMP SHARED(rhtdz,rhtdzz,rxydz,rxydzz,l_noco,dzmagmom,ddzmagmom,idx) &
       !$OMP SHARED(ifftd2,rho,grad) &
       !$OMP PRIVATE(ip,js,iq,cqpw,bf2,rhti,rhdx,rhdy,rhdz,rhdxx,rhdyy,rhdzz) &
       !$OMP PRIVATE(rhdxy,rhdzx,rhdyz,dxmagmom,dymagmom,ddxmagmom,ddymagmom) &
       !$OMP PRIVATE(chdens,idx_loc)
       ALLOCATE ( rhdx(0:ifftd2-1,jspins),rhdy(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdz(0:ifftd2-1,jspins),rhdxx(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdyy(0:ifftd2-1,jspins),rhdzz(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdyz(0:ifftd2-1,jspins),rhdzx(0:ifftd2-1,jspins) )
       ALLOCATE ( rhdxy(0:ifftd2-1,jspins))
       ALLOCATE ( cqpw(stars%ng2-1))
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

                DO iq=1,stars%ng2-1
                   cqpw(iq)=ImagUnit*vacxy(ip,iq,ivac,js)
                ENDDO

                rhti = 0.0                    ! d(rho)/atoms%dx is obtained by a FFT of i*gx*vacxy
                ! (vacz is set to zero and gx is included in
                !    dn/atoms =  FFT(0,i*gx*vacxy)



                CALL fft2d(stars, rhdx(0,js),bf2, zro,rhti,cqpw,+1,firstderiv=[1.,0.,0.],cell=cell)
                !TODO    &                 pgft2x)

                rhti = 0.0
                CALL fft2d(    &               ! dn/dy =  FFT(0,i*gy*vacxy)&
                      stars, rhdy(0,js),bf2, zro,rhti,cqpw, +1,firstderiv=[0.,1.,0.],cell=cell)

                rhti = 0.0
                CALL fft2d(     &              ! dn/dz = FFT(rhtdz,rxydz)&
                        stars, rhdz(0,js),bf2, rhtdz(ip,js),rhti,rxydz(ip,:,js), +1)

                DO iq=1,stars%ng2-1
                   cqpw(iq)=-vacxy(ip,iq,ivac,js)
                ENDDO

                rhti = 0.0
                CALL fft2d(      &          ! d2n/dx2 = FFT(0,-gx^2*vacxy)&
                       stars, rhdxx(0,js),bf2, zro,rhti,cqpw, +1,firstderiv=[1.0,0.,0.],secondderiv=[1.0,0.,0.],cell=cell)

                rhti = 0.0
                CALL fft2d(       &          ! d2n/dy2 = FFT(0,-gy^2*vacxy)&
                      stars, rhdyy(0,js),bf2, zro,rhti,cqpw, +1,firstderiv=[0.,1.0,0.],secondderiv=[0.,1.0,0.],cell=cell)

                rhti = 0.0
                CALL fft2d(        &         ! d2n/dz2 = FFT(rhtdzz,rxydzz)&
                       stars, rhdzz(0,js),bf2, rhtdzz(ip,js),rhti,rxydzz(ip,:,js), +1)


                DO iq=1,stars%ng2-1
                   cqpw(iq)=ImagUnit*rxydz(ip,iq,js)
                ENDDO

                rhti = 0.0
                CALL fft2d(         &         ! d2n/dyz = FFT(0,i*gy*rxydz)&
                       stars, rhdyz(0,js),bf2, zro,rhti,cqpw, +1,firstderiv=[0.,1.0,0.],cell=cell)

                rhti = 0.0
                CALL fft2d(          &        ! d2n/dzx = FFT(0,i*gx*rxydz)&
                       stars, rhdzx(0,js),bf2, zro,rhti,cqpw, +1,firstderiv=[1.,0.0,0.],cell=cell)

                DO iq=1,stars%ng2-1
                   cqpw(iq)=-vacxy(ip,iq,ivac,js)
                ENDDO

                rhti = 0.0
                CALL fft2d(           &    ! d2n/dxy = FFT(0,-gx*gy*vacxy)&
                      stars, rhdxy(0,js),bf2, zro,rhti,cqpw, +1,firstderiv=[0.,1.0,0.],secondderiv=[1.,0.0,0.],cell=cell)

             END DO ! js=1,jspins


             IF (l_noco) THEN
                ! ! In non-collinear calculations the derivatives of |m| are calculated
                ! ! in real-space. The derivatives of the charge density, that are
                ! ! already calculated in g-space, will be used.

                CALL grdrsvac(magmom(0,ip),cell%bmat,3*stars%mx1,3*stars%mx2,fixed_ndvgrd, dxmagmom,dymagmom)
                DO i=0,9*stars%mx1*stars%mx2-1
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
                DO i=0,9*stars%mx1*stars%mx2-1
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
                DO i=0,9*stars%mx1*stars%mx2-1
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
!          if(oneD%odi%d1)then
!          rd = cell%z1 + vacuum%delz*(ip-1)
!!$             CALL od_mkgxyz3(&
!!$                  &           ifftd2,input%jspins,ifftd2,input%jspins,&
!!$                  &           af2,rd,rhdx,rhdy,rhdz,rhdxx,rhdyy,&
!!$                  &           rhdzz,rhdyz,rhdzx,rhdxy,&
!!$                  &           agr,agru,agrd,g2r,g2ru,g2rd,&
!!$                  &           gggr,gggru,gggrd,gzgr)
!             CALL judft_error("OneD not implemented")
!          ELSE
             idx_loc = idx + (ip-1)* ifftd2
             CALL mkgxyz3(rho(idx_loc:idx_loc+ifftd2-1,:),rhdx,rhdy, rhdz,rhdxx,rhdyy,rhdzz,rhdyz,rhdzx,rhdxy, idx_loc-1,grad)
!          endif

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
                rho(idx+ip-nmz0,js)= vacz(ip,ivac,js)
             END DO
          ELSE
             mx(0) = vacz(ip,ivac,3)
             my(0) = vacz(ip,ivac,4)
             chdens= (vacz(ip,ivac,1)+vacz(ip,ivac,2))/2.
             magmom(0,1)= mx(0)**2 + my(0)**2 + ((vacz(ip,ivac,1)-vacz(ip,ivac,2))/2.)**2
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


!          if(oneD%odi%d1)then
!             CALL od_mkgz(&
!                              cell%z1,vacuum%nmzxy,vacuum%delz,&
!                              nmzdiff,jspins,&
!                              rhtz(vacuum%nmzxy+1,1),rhtz(vacuum%nmzxy+1,jspins),&
!                              rhtdz(vacuum%nmzxy+1,1), rhtdz(vacuum%nmzxy+1,jspins),&
!                              rhtdzz(vacuum%nmzxy+1,1),rhtdzz(vacuum%nmzxy+1,jspins),&
!                              agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,&
!                              gzgr)
!             CALL judft_error("OneD not implemented")
!          ELSE
             CALL mkgz(nmzdiff,jspins, rho(nmz0:,1),rho(nmz0:,jspins),&
             rhtdz(nmz0:,1),rhtdz(nmz0:,jspins),rhtdzz(nmz0:,1),&
                  rhtdzz(nmz0:,jspins),idx,grad)

!          endif
       ENDIF

       !       calculate vxc for z now beyond warping region
       idx=idx+nmzdiff
    ENDDO    ! loop over vacua (ivac)



  END SUBROUTINE vac_to_grid

  subroutine vac_from_grid(stars,vacuum,v_xc,ifft2d,vacz,vacxy)
    use m_types_stars
    use m_types_vacuum
    use m_fft2d
    type(t_stars),intent(in)  :: stars
    type(t_vacuum),intent(in) :: vacuum

    real, INTENT(IN)          :: v_xc(:,:)
    INTEGER,INTENT(IN)        :: ifft2d
    real, intent(INOUT)       :: vacz(:,:,:)
    complex,intent(INOUT)     :: vacxy(:,:,:,:)

    REAL                    :: fgz,rhti
    COMPLEX, ALLOCATABLE    :: fgxy(:)
    REAL, ALLOCATABLE       :: bf2(:)
    INTEGER                 :: js,irec2,idx,ivac,ip

    ALLOCATE ( fgxy(stars%ng2-1),bf2(ifft2d) )

    DO js = 1,size(v_xc,2)
      idx=1
      DO ivac=1,vacuum%nvac
        DO ip=1,vacuum%nmzxy
          !
          !           ----> 2-d back fft to g space
          !
          bf2=0.0
          CALL fft2d(stars, v_xc(idx:idx-1+ifft2d,js),bf2, fgz,rhti,fgxy, -1)
          idx=idx+ifft2d
          !            ----> and add vxc to coulomb potential
          !                  the g||.eq.zero component is added to vxc%vacz
          !
          !vxc%vacz(ip,ivac,js) = fgz + vxc%vacz(ip,ivac,js)
          vacz(ip,ivac,js) = fgz + vacz(ip,ivac,js)
          !
          !            the g||.ne.zero components are added to vxc%vacxy
          !
          DO irec2 = 1,stars%ng2-1
            vacxy(ip,irec2,ivac,js)=vacxy(ip,irec2,ivac,js)+fgxy(irec2)
          ENDDO
        enddo

        nmz0= max(1,vacuum%nmzxy+1+(fixed_ndvgrd/2)-fixed_ndvgrd)

        DO ip = nmz0,vacuum%nmz
          if (ip>vacuum%nmzxy)vacz(ip,ivac,js) = vacz(ip,ivac,js) + v_xc(idx,js)
          idx=idx+1
        ENDDO
      END DO
    ENDDO

  end subroutine
end module
