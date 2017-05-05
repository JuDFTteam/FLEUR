MODULE m_vvacxcg
  use m_juDFT
  !-----------------------------------------------------------------------
  !     calculates 2-d star function coefficients of exchange-correlation*
  !     potential in the vacuum regions and adds them to the corresponding
  !     coeffs of the coulomb potential            c.l.fu, r.podloucky   *
  !     for the gradient contribution.   t.a. 1996
  !-----------------------------------------------------------------------
CONTAINS
  SUBROUTINE vvacxcg(&
       &           ifftd2,stars,vacuum,noco,oneD,&
       &           cell,xcpot,input,obsolete,&
       &           ichsmrg,&
       &           rhtxy,rht,cdomvxy,cdomvz,&
       &           vxy,vz,rhmn,&
       &           excxy,excz)

    !-----------------------------------------------------------------------
    !     instead of vvacxcor.f: the different exchange-correlation
    !     potentials defined through the key icorr are called through
    !     the driver subroutine vxcallg.f, subroutines vectorized
    !     in case of total = .true. calculates the ex-corr. energy
    !     density through the driver subroutine excallg.f
    !     ** r.pentcheva 08.05.96
    !-----------------------------------------------------------------------
    USE m_grdrsvac
    USE m_grdchlh
    USE m_mkgz
    USE m_mkgxyz3
    USE m_od_mkgxyz3
    USE m_od_mkgz
    USE m_fft2d
    USE m_xcallg, ONLY : vxcallg,excallg
    USE m_types
    IMPLICIT NONE
    TYPE(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_obsolete),INTENT(IN)  :: obsolete
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ifftd2
    INTEGER, INTENT (INOUT) :: ichsmrg
    REAL,    INTENT (INOUT) :: rhmn
    !-odim
    !+odim
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rht(vacuum%nmzd,2,input%jspins)
    COMPLEX, INTENT (IN) :: rhtxy(vacuum%nmzxy,stars%ng2-1,2,input%jspins)
    COMPLEX, INTENT (IN) :: cdomvz(vacuum%nmzd,2)
    COMPLEX, INTENT (IN) :: cdomvxy(vacuum%nmzxy,stars%ng2-1,2)
    REAL,    INTENT (OUT) :: excz(vacuum%nmzd,2)
    COMPLEX, INTENT (OUT) :: excxy(vacuum%nmzxy,stars%ng2-1,2)
    REAL,    INTENT (INOUT) :: vz(vacuum%nmzd,2,input%jspins)
    COMPLEX, INTENT (INOUT) :: vxy(vacuum%nmzxy,stars%ng2-1,2,input%jspins)
    !     ..
    !     .. Local Scalars ..
    INTEGER :: js,nt,i,iq,irec2,nmz0,nmzdiff,ivac,ip
    REAL    :: rhti,zro,fgz,rhmnv,d_15,bmat1(3,3),rd
    COMPLEX :: ci 
    LOGICAL :: lwbc              ! if true, white-bird trick.
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: af2(:,:),bf2(:),agr(:),agru(:),agrd(:),g2r(:)
    REAL, ALLOCATABLE :: g2ru(:),g2rd(:),gggr(:),gggru(:),gggrd(:)
    REAL, ALLOCATABLE :: gzgr(:),rhdx(:,:),rhdy(:,:),rhdz(:,:)
    REAL, ALLOCATABLE :: rhdxx(:,:),rhdyy(:,:),rhtdz(:,:),rhtdzz(:,:)
    REAL, ALLOCATABLE :: rhdzz(:,:),rhdyz(:,:),rhdzx(:,:),rhdxy(:,:)
    REAL, ALLOCATABLE :: vx(:,:),vxc(:,:),exc(:),vxz(:,:),vxcz(:,:)
    REAL, ALLOCATABLE :: rxydzr(:),rxydzi(:)
    REAL, ALLOCATABLE :: rxydzzr(:),rxydzzi(:),rhtxyr(:),rhtxyi(:)
    REAL, ALLOCATABLE :: rhtxc(:,:),rhtz(:,:),dummy(:)
    COMPLEX, ALLOCATABLE :: fgxy(:),rxydz(:,:,:),rxydzz(:,:,:),cqpw(:)
    !     ..
    !     for the noco-case only 
    REAL :: chdens
    REAL, ALLOCATABLE :: magmom(:,:),&
         &                     dxmagmom(:),ddxmagmom(:,:),&
         &                     dymagmom(:),ddymagmom(:,:), &
         &                     dzmagmom(:,:),ddzmagmom(:,:)
    REAL, ALLOCATABLE :: mx(:),my(:) 
    REAL, ALLOCATABLE :: af2noco(:,:,:)

    !     .. unused input (needed for other noco GGA-implementations) ..



    lwbc     = .false.
    d_15     = 1.e-15
    zro      = 0.0
    ci       = cmplx(0.,1.)
    nt       = ifftd2
    if(oneD%odi%d1)then
       bmat1(:,:) = 0.
       bmat1(1,1) = cell%bmat(3,3)
       bmat1(2,2) = 1.
    else
       bmat1(:,:) = cell%bmat(:,:)
    endif

    WRITE (6,'(/'' ifftd2,vacuum%nmz='',2i7)') ifftd2,vacuum%nmz
    WRITE(6,'('' 9990nmzxy='',2i5)') vacuum%nmzxy

    ALLOCATE ( rxydz(vacuum%nmzxy,stars%ng2-1,input%jspins),rxydzz(vacuum%nmzxyd,stars%ng2-1,input%jspins) )
    ALLOCATE ( rhtdz(vacuum%nmzd,input%jspins),rhtdzz(vacuum%nmzd,input%jspins) )

    DO ivac=1,vacuum%nvac 

       ! the charge density in vacuum is expanded in 2-dim stars on a mesh
       ! in z-direction. the g||.ne.zero-components expand from 1 to nmzxy
       ! the g||.eq.zero-components expand from 1 to nmz
       ! first we calculate vxc in the warping region 

       IF (noco%l_noco) THEN

          ALLOCATE ( magmom(0:ifftd2-1,vacuum%nmzxy) ) 
          ALLOCATE ( dzmagmom(0:ifftd2-1,vacuum%nmzxy) ) 
          ALLOCATE ( ddzmagmom(0:ifftd2-1,vacuum%nmzxy) ) 
          ALLOCATE ( mx(0:ifftd2-1),my(0:ifftd2-1) )
          ALLOCATE ( dxmagmom(0:ifftd2-1),dymagmom(0:ifftd2-1) )
          ALLOCATE ( ddxmagmom(0:ifftd2-1,2),ddymagmom(0:ifftd2-1,2) ) 
          ALLOCATE ( af2noco(0:ifftd2-1,vacuum%nmzxy,input%jspins),bf2(0:ifftd2-1) ) 

          ! Transform charge and magnetization to real-space.
          ! In the collinear case that is done later within
          ! another loop over the vacuum-layers in order to 
          ! save memory.

          DO ip=1,vacuum%nmzxy

             DO js=1,input%jspins
                CALL fft2d(&
                     &                   stars,&
                     &                   af2noco(0,ip,js),bf2,&
                     &                   rht(ip,ivac,js),0.,rhtxy(ip,1,ivac,js),&
                     &                   vacuum%nmzxy,+1)
             END DO
             CALL fft2d(&
                  &                 stars,&
                  &                 mx,my,&
                  &                 REAL(cdomvz(ip,ivac)),AIMAG(cdomvz(ip,ivac)),&
                  &                 cdomvxy(ip,1,ivac),&
                  &                 vacuum%nmzxy,+1)

             DO i=0,9*stars%mx1*stars%mx2-1
                magmom(i,ip)= mx(i)**2 + my(i)**2 +&
                     &          ((af2noco(i,ip,1)-af2noco(i,ip,2))/2.)**2
                magmom(i,ip)= SQRT(magmom(i,ip))
                chdens= af2noco(i,ip,1)/2.+af2noco(i,ip,2)/2.
                af2noco(i,ip,1)= chdens + magmom(i,ip)
                af2noco(i,ip,2)= chdens - magmom(i,ip)
             END DO

          END DO ! ip=1,vacuum%nmzxy 
          DEALLOCATE ( bf2 )
       END IF ! noco%l_noco 

       !      ENDDO    ! ivac
       !      DO ivac = 1,nvac

       IF (xcpot%igrd.GT.0) THEN
          DO js=1,input%jspins
             !
             ! calculate first (rhtdz) & second (rhtdzz) derivative of rht(1:nmz)
             !
             ALLOCATE ( dummy(vacuum%nmz) )
             CALL grdchlh(&
                  &                   0,1,vacuum%nmz,vacuum%delz,dummy,rht(1,ivac,js),obsolete%ndvgrd,&
                  &                   rhtdz(1,js),rhtdzz(1,js))
             DEALLOCATE ( dummy )
             ALLOCATE ( rhtxyr(vacuum%nmzxy), rhtxyi(vacuum%nmzxy),dummy(vacuum%nmzxy) )
             ALLOCATE ( rxydzr(vacuum%nmzxy), rxydzi(vacuum%nmzxy) )
             ALLOCATE ( rxydzzr(vacuum%nmzxy),rxydzzi(vacuum%nmzxy) )

             DO iq = 1, stars%ng2-1
                !
                ! calculate first (rxydz) & second (rxydzz) derivative of rhtxy:
                !
                DO ip=1,vacuum%nmzxy
                   rhtxyr(ip)=rhtxy(ip,iq,ivac,js)
                ENDDO
                CALL grdchlh(&
                     &                     0,1,vacuum%nmzxy,vacuum%delz,dummy,rhtxyr,obsolete%ndvgrd,&
                     &                     rxydzr,rxydzzr) 

                DO ip=1,vacuum%nmzxy
                   rhtxyi(ip)=aimag(rhtxy(ip,iq,ivac,js))
                ENDDO
                CALL grdchlh(&
                     &                     0,1,vacuum%nmzxy,vacuum%delz,dummy,rhtxyi,obsolete%ndvgrd,&
                     &                     rxydzi,rxydzzi)

                DO ip=1,vacuum%nmzxy
                   rxydz(ip,iq,js)=cmplx(rxydzr(ip),rxydzi(ip))
                   rxydzz(ip,iq,js)=cmplx(rxydzzr(ip),rxydzzi(ip))
                ENDDO

             ENDDO ! loop over 2D stars (iq)

             DEALLOCATE ( rhtxyr,rhtxyi,rxydzr,rxydzi,rxydzzr,rxydzzi )
             DEALLOCATE ( dummy )

          ENDDO ! input%jspins

          IF (noco%l_noco) THEN 
             !  calculate  dzmagmom = d magmom / d z  and ddzmagmom= d dmagmom / d z 

             ALLOCATE ( rhtxyr(vacuum%nmzxy),dummy(vacuum%nmzxy)   )
             ALLOCATE ( rxydzr(vacuum%nmzxy),rxydzzr(vacuum%nmzxy) )
             DO i=0,9*stars%mx1*stars%mx2-1 
                DO ip=1,vacuum%nmzxy
                   rhtxyr(ip)=magmom(i,ip)
                ENDDO
                CALL grdchlh(&
                     &                     0,1,vacuum%nmzxy,vacuum%delz,dummy,rhtxyr,obsolete%ndvgrd,&
                     &                     rxydzr,rxydzzr)
                DO ip=1,vacuum%nmzxy
                   dzmagmom(i,ip)= rxydzr(ip)
                   ddzmagmom(i,ip)= rxydzzr(ip)
                ENDDO
             END DO
             DEALLOCATE ( rhtxyr,rxydzr,rxydzzr,dummy )
          END IF ! noco%l_noco 

       ENDIF   ! xcpot%igrd.GT.0

       !       WRITE(6,'('' 9990nmzxy='',2i5)') nmzxy
       ALLOCATE ( rhdx(0:ifftd2-1,input%jspins),rhdy(0:ifftd2-1,input%jspins) )
       ALLOCATE ( rhdz(0:ifftd2-1,input%jspins),rhdxx(0:ifftd2-1,input%jspins) )
       ALLOCATE ( rhdyy(0:ifftd2-1,input%jspins),rhdzz(0:ifftd2-1,input%jspins) )
       ALLOCATE ( rhdyz(0:ifftd2-1,input%jspins),rhdzx(0:ifftd2-1,input%jspins) )
       ALLOCATE ( rhdxy(0:ifftd2-1,input%jspins),gggrd(0:ifftd2-1) )
       ALLOCATE ( agr(0:ifftd2-1),agru(0:ifftd2-1),agrd(0:ifftd2-1) )
       ALLOCATE ( g2r(0:ifftd2-1),g2ru(0:ifftd2-1),g2rd(0:ifftd2-1) )
       ALLOCATE ( gggr(0:ifftd2-1),gggru(0:ifftd2-1),gzgr(0:ifftd2-1) )
       ALLOCATE ( vxc(0:ifftd2-1,input%jspins),exc(0:ifftd2-1) )
       ALLOCATE ( cqpw(stars%ng2-1),af2(0:ifftd2-1,input%jspins) )
       ALLOCATE ( fgxy(stars%ng2-1),bf2(0:ifftd2-1) )

       ALLOCATE( vx(0:ifftd2-1,input%jspins) )

 
       rd = 0.0
       af2=0.0
       DO ip = 1,vacuum%nmzxy
          ! loop over warping region

          IF (.not. noco%l_noco) THEN
             ! Transform charge and magnetization to real-space.

             DO js=1,input%jspins
                CALL fft2d(&
                     &               stars,&
                     &               af2(0,js),bf2,&
                     &               rht(ip,ivac,js),0.,rhtxy(ip,1,ivac,js),&
                     &               vacuum%nmzxyd,+1)
             END DO

          ELSE

             DO i=0,9*stars%mx1*stars%mx2-1
                af2(i,1)= af2noco(i,ip,1) 
                af2(i,2)= af2noco(i,ip,2)
             END DO

          END IF

          IF (xcpot%igrd > 0) THEN 
             ! calculate derivatives with respect to x,y in g-space 
             ! and transform them to real-space.  

             DO js = 1,input%jspins

                DO iq=1,stars%ng2-1
                   cqpw(iq)=ci*rhtxy(ip,iq,ivac,js)
                ENDDO

                rhti = 0.0                    ! d(rho)/atoms%dx is obtained by a FFT of i*gx*rhtxy
                ! (rht is set to zero and gx is included in 
                !    dn/atoms =  FFT(0,i*gx*rhtxy)

              

                CALL fft2d(          &
                     &stars,        &
                     &rhdx(0,js),bf2,&
                     &zro,rhti,cqpw,       &
                     &1,+1,stars%ft2_gfx)
                !TODO    &                 pgft2x) 

                rhti = 0.0
                CALL fft2d(    &               ! dn/dy =  FFT(0,i*gy*rhtxy)&
                     & stars,&
                     & rhdy(0,js),bf2,&
                     & zro,rhti,cqpw,&
                     & 1,+1,stars%ft2_gfy)

                rhti = 0.0
                CALL fft2d(     &              ! dn/dz = FFT(rhtdz,rxydz)&
                     &   stars,&
                     &   rhdz(0,js),bf2,&
                     &   rhtdz(ip,js),rhti,rxydz(ip,1,js),&
                     &   vacuum%nmzxyd,+1)


                DO iq=1,stars%ng2-1
                   cqpw(iq)=-rhtxy(ip,iq,ivac,js)
                ENDDO

                rhti = 0.0
                CALL fft2d(      &          ! d2n/dx2 = FFT(0,-gx^2*rhtxy)&
                     &  stars,&
                     &  rhdxx(0,js),bf2,&
                     &  zro,rhti,cqpw,&
                     &  1,+1,stars%ft2_gfx*stars%ft2_gfx)

                rhti = 0.0
                CALL fft2d(       &          ! d2n/dy2 = FFT(0,-gy^2*rhtxy)&
                     & stars,&
                     & rhdyy(0,js),bf2,&
                     & zro,rhti,cqpw,&
                     & 1,+1,stars%ft2_gfy*stars%ft2_gfy)

                rhti = 0.0
                CALL fft2d(        &         ! d2n/dz2 = FFT(rhtdzz,rxydzz)&
                     &  stars,&
                     &  rhdzz(0,js),bf2,&
                     &  rhtdzz(ip,js),rhti,rxydzz(ip,1,js),&
                     &  vacuum%nmzxyd,+1)


                DO iq=1,stars%ng2-1
                   cqpw(iq)=ci*rxydz(ip,iq,js)
                ENDDO

                rhti = 0.0
                CALL fft2d(         &         ! d2n/dyz = FFT(0,i*gy*rxydz)&
                     &  stars,&
                     &  rhdyz(0,js),bf2,&
                     &  zro,rhti,cqpw,&
                     &  1,+1,stars%ft2_gfy)

                rhti = 0.0
                CALL fft2d(          &        ! d2n/dzx = FFT(0,i*gx*rxydz)&
                     &  stars,&
                     &  rhdzx(0,js),bf2,&
                     &  zro,rhti,cqpw,&
                     &  1,+1,stars%ft2_gfx)

                DO iq=1,stars%ng2-1
                   cqpw(iq)=-rhtxy(ip,iq,ivac,js)
                ENDDO

                rhti = 0.0
                CALL fft2d(           &    ! d2n/dxy = FFT(0,-gx*gy*rhtxy)&
                     & stars,&
                     & rhdxy(0,js),bf2,&
                     & zro,rhti,cqpw,&
                     & 1,+1,stars%ft2_gfy*stars%ft2_gfx)

             END DO ! js=1,input%jspins


             IF (noco%l_noco) THEN
                ! ! In non-collinear calculations the derivatives of |m| are calculated
                ! ! in real-space. The derivatives of the charge density, that are 
                ! ! already calculated in g-space, will be used. 

                CALL grdrsvac(&
                     &               magmom(0,ip),bmat1,3*stars%mx1,3*stars%mx2,obsolete%ndvgrd,&
                     &               dxmagmom,dymagmom) 
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

                CALL grdrsvac(&
                     &               dxmagmom,bmat1,3*stars%mx1,3*stars%mx2,obsolete%ndvgrd, &
                     &               ddxmagmom(0,1),ddymagmom(0,1))
                CALL grdrsvac(&
                     &               dymagmom,bmat1,3*stars%mx1,3*stars%mx2,obsolete%ndvgrd, &
                     &               ddxmagmom(0,2),ddymagmom(0,2))
                DO i=0,9*stars%mx1*stars%mx2-1
                   chdens= rhdxx(i,1)/2.+rhdxx(i,2)/2. 
                   rhdxx(i,1)= chdens + ddxmagmom(i,1) 
                   rhdxx(i,2)= chdens - ddxmagmom(i,1) 
                   chdens= rhdyy(i,1)/2.+rhdyy(i,2)/2. 
                   rhdyy(i,1)= chdens + ddymagmom(i,2) 
                   rhdyy(i,2)= chdens - ddymagmom(i,2) 
                   chdens= rhdxy(i,1)/2.+rhdxy(i,2)/2. 
                   rhdxy(i,1)= chdens + &
                        &           ( ddxmagmom(i,2) + ddymagmom(i,1) )/2.
                   rhdxy(i,2)= chdens - &
                        &           ( ddxmagmom(i,2) + ddymagmom(i,1) )/2.
                END DO
                CALL grdrsvac(&
                     &               dzmagmom(0,ip),bmat1,3*stars%mx1,3*stars%mx2,obsolete%ndvgrd, &
                     &               ddxmagmom(0,1),ddymagmom(0,1))
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

             END IF ! noco%l_noco  

          END IF ! xcpot%igrd > 0 
          !
          ! set minimal value of af2 to 1.0e-13
          !

          af2=max(af2,10e-13)


          !   calculate the quantities such as abs(grad(rho)),.. used in
          !c  evaluating the gradient contributions to potential and energy.

          rd = cell%z1 + vacuum%delz*(ip-1)

          if(oneD%odi%d1)then
             CALL od_mkgxyz3(&
                  &           xcpot%igrd,ifftd2,input%jspins,ifftd2,input%jspins,&
                  &           af2,rd,rhdx,rhdy,rhdz,rhdxx,rhdyy,&
                  &           rhdzz,rhdyz,rhdzx,rhdxy,&
                  &           agr,agru,agrd,g2r,g2ru,g2rd,&
                  &           gggr,gggru,gggrd,gzgr)
          else
             CALL mkgxyz3(&
                  &           xcpot%igrd,ifftd2,input%jspins,ifftd2,input%jspins,af2,rhdx,rhdy,&
                  &           rhdz,rhdxx,rhdyy,rhdzz,rhdyz,rhdzx,rhdxy,&
                  &           agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,&
                  &           gzgr)
          endif

          !     rhmnv: rho_minimum_vacuum.

          rhmnv=10.e+10
          DO js=1,input%jspins
             DO i=0,stars%kimax2
                af2(i,js)=max(af2(i,js),d_15)
                rhmnv=min(rhmnv,af2(i,js))
             ENDDO
          ENDDO

          IF (rhmnv.lt.rhmn) THEN
             rhmn=rhmnv
             ichsmrg=3
          ENDIF

          IF (rhmn.LT.obsolete%chng) THEN
             WRITE(6,'(/'' rhmn.lt.obsolete%chng. rhmn,obsolete%chng='',2d9.2)') rhmn,obsolete%chng
             !             CALL juDFT_error("vvacxcg: rhmn.lt.chng",calledby="vvacxcg")
          ENDIF

          !         calculate the exchange-correlation potential in  real space
          !

          CALL vxcallg(&
               &                 xcpot%icorr,lwbc,input%jspins,nt,nt,af2,agr,agru,agrd,&
               &                 g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,&
               &                 vx,vxc)



          DO js = 1,input%jspins
             !
             !           ----> 2-d back fft to g space
             !
             bf2=0.0
             CALL fft2d(&
                  &                 stars,&
                  &                 vxc(0,js),bf2,&
                  &                 fgz,rhti,fgxy,&
                  &                 1,-1)

             !            ----> and add vxc to coulomb potential
             !                  the g||.eq.zero component is added to vz
             !
             vz(ip,ivac,js) = fgz + vz(ip,ivac,js)
             !
             !            the g||.ne.zero components are added to vxy
             !
             DO irec2 = 1,stars%ng2-1
                vxy(ip,irec2,ivac,js)=vxy(ip,irec2,ivac,js)+fgxy(irec2)
             ENDDO

          END DO


          !         calculate the exchange-correlation energy density in  real space
          !
          IF (input%total) THEN

             CALL excallg(&
                  &                   xcpot%icorr,lwbc,input%jspins,nt,af2,agr,agru,agrd,&
                  &                   g2r,g2ru,g2rd, gggr,gggru,gggrd,gzgr,&
                  &                   exc)

             !           ----> 2-d back fft to g space
             !
             bf2=0.0
             CALL fft2d(&
                  &                 stars,&
                  &                 exc,bf2,&
                  &                 excz(ip,ivac),rhti,excxy(ip,1,ivac),&
                  &                 vacuum%nmzxyd,-1)

          ENDIF

       END DO ! ip=1,vacuum%nmzxy 
       DEALLOCATE ( rhdx,rhdy,rhdz,rhdxx,rhdyy,rhdzz )
       DEALLOCATE ( cqpw,fgxy,     rhdyz,rhdzx,rhdxy )

       IF (noco%l_noco) THEN
          DEALLOCATE ( dzmagmom,ddzmagmom,dxmagmom,af2noco )
          DEALLOCATE ( dymagmom,ddxmagmom,ddymagmom )
       END IF

       ! now treat the non-warping region 

       nmzdiff = vacuum%nmz - vacuum%nmzxy
       !       WRITE(6,'(/'' 9992excz''/(8f15.7))') (excz(ip,1),ip=1,nmz)
       WRITE(6,'(/'' 9992nmzdiff='',i5)') nmzdiff

       ! The non-warping region runs from nmzxy+1 to nmz.
       ! The values from nmz0 to nmzxy are taken into account in order
       ! to get the real-space derivative smooth around nmzxy+1. 
       nmz0= vacuum%nmzxy+1+(obsolete%ndvgrd/2)-obsolete%ndvgrd
       IF (nmz0 <= 0) THEN ! usually vacuum%nmzxy>obsolete%ndvgrd 
          nmz0= 1
       END IF

       ALLOCATE ( rhtz(vacuum%nmzd,input%jspins) )

       DO ip=nmz0,vacuum%nmz 
          IF (.not. noco%l_noco) THEN
             DO js=1,input%jspins 
                rhtz(ip,js)= rht(ip,ivac,js)
             END DO
          ELSE
             af2(0,1) = rht(ip,ivac,1)
             af2(0,2) = rht(ip,ivac,2)
             mx(0)= REAL(cdomvz(ip,ivac))
             my(0)= AIMAG(cdomvz(ip,ivac))
             chdens= (af2(0,1)+af2(0,2))/2.
             magmom(0,1)= mx(0)**2 + my(0)**2 +&
                  &                   ((af2(0,1)-af2(0,2))/2.)**2
             magmom(0,1)= SQRT(magmom(0,1))
             rhtz(ip,1)= chdens + magmom(0,1)
             rhtz(ip,2)= chdens - magmom(0,1) 
          END IF
       END DO

       IF (noco%l_noco) THEN 
          DEALLOCATE ( magmom,mx,my ) 
          ALLOCATE ( dummy(vacuum%nmz) )
          DO js=1,input%jspins
             CALL grdchlh(&
                  &                   0,1,vacuum%nmz-nmz0+1,vacuum%delz,dummy,rhtz(nmz0,js),obsolete%ndvgrd,&
                  &                   rhtdz(nmz0,js),rhtdzz(nmz0,js))
          END DO
          DEALLOCATE ( dummy )

       END IF

       !       calculate the quantities such as abs(grad(rho)),.. used in
       !c      evaluating the gradient contributions to potential and
       !c      energy.

       agr(:)=0.0 ; agru(:)=0.0 ; agrd(:)=0.0 ; g2r(:)=0.0
       g2ru(:)=0.0 ; g2rd(:)=0.0 ; gggr(:)=0.0 ; gggru(:)=0.0
       gggrd(:)=0.0 ; gzgr(:)=0.0

       IF (xcpot%igrd.gt.0)  THEN
          if(oneD%odi%d1)then
             CALL od_mkgz(&
                  &            cell%z1,vacuum%nmzxy,vacuum%delz,&
                  &            nmzdiff,input%jspins,&
                  &            rhtz(vacuum%nmzxy+1,1),rhtz(vacuum%nmzxy+1,input%jspins),&
                  &            rhtdz(vacuum%nmzxy+1,1), rhtdz(vacuum%nmzxy+1,input%jspins),&
                  &            rhtdzz(vacuum%nmzxy+1,1),rhtdzz(vacuum%nmzxy+1,input%jspins),&
                  &            agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,&
                  &            gzgr)

          else
             CALL mkgz(nmzdiff,input%jspins,&
                  &            rhtz(vacuum%nmzxy+1,1),rhtz(vacuum%nmzxy+1,input%jspins),&
                  &            rhtdz(vacuum%nmzxy+1,1),rhtdz(vacuum%nmzxy+1,input%jspins),&
                  &            rhtdzz(vacuum%nmzxy+1,1),rhtdzz(vacuum%nmzxy+1,input%jspins),&
                  &            agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,&
                  &            gzgr)
          endif
       ENDIF

       !       calculate vxc for z now beyond warping region
       DEALLOCATE ( af2)
       ALLOCATE ( rhtxc(vacuum%nmzd,input%jspins),vxcz(vacuum%nmzd,input%jspins) )
       ALLOCATE ( vxz(vacuum%nmzd,input%jspins) )

       DO js=1,input%jspins
          !DO i=0,ifftd2-1
          !  af2(i,js)=0.0
          !ENDDO
          DO ip=vacuum%nmzxy+1,vacuum%nmz
             rhtxc(ip-vacuum%nmzxy,js) = max(rhtz(ip,js),d_15) !+gb
             rhmnv=min(rhmnv,rhtxc(ip-vacuum%nmzxy,js))
          ENDDO
       ENDDO

       IF(rhmnv.lt.rhmn) THEN
          rhmn=rhmnv
          ichsmrg=4
       ENDIF

       IF (rhmn.lt.obsolete%chng) THEN
          WRITE (6,'('' rhmn.lt.obsolete%chng. rhmn,obsolete%chng='',2d9.2)') rhmn,obsolete%chng
          !           CALL juDFT_error("vvacxcg: rhmn.lt.chng",calledby="vvacxcg")
       ENDIF

       CALL vxcallg(&
            &             xcpot%icorr,lwbc,input%jspins,vacuum%nmzd,nmzdiff,rhtxc,agr(:vacuum%nmzd),agru,&
            &               agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,&
            &               vxz,vxcz)

       DO js = 1,input%jspins
          DO ip = vacuum%nmzxy + 1,vacuum%nmz
             vz(ip,ivac,js) = vz(ip,ivac,js) + vxcz(ip-vacuum%nmzxy,js)
          ENDDO
       ENDDO

       !
       WRITE (6,fmt=8020) ivac, (vz(vacuum%nmz,ivac,js),js=1,input%jspins)
       WRITE(16,fmt=8020) ivac, (vz(vacuum%nmz,ivac,js),js=1,input%jspins)
8020   FORMAT(/,5x,'vacuum zero for vacuum',i3,' = ',2f14.10)
       !
       !     calculate the ex-corr. energy density now beyond warping region
       !
       IF (input%total) THEN
          CALL excallg(&
               &                   xcpot%icorr,lwbc,input%jspins,nmzdiff,rhtxc,agr,agru,&
               &                   agrd,g2r,g2ru,g2rd,gggr,gggru,gggrd,gzgr,&
               &                   excz(vacuum%nmzxy+1,ivac))
       ENDIF


       DEALLOCATE ( bf2)
       DEALLOCATE ( agr,agru )
       DEALLOCATE ( agrd,g2r )
       DEALLOCATE ( g2ru,g2rd )
       DEALLOCATE ( gggr,gggru )
       DEALLOCATE ( gggrd,gzgr,vx,vxc,exc,vxz,vxcz,rhtz,rhtxc )

    ENDDO    ! loop over vacua (ivac)
    DEALLOCATE ( rhtdz,rhtdzz,rxydz,rxydzz )



  END SUBROUTINE vvacxcg
END MODULE m_vvacxcg
