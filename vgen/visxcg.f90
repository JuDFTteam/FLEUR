!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_visxcg
  USE m_juDFT
  !     ******************************************************
  !     subroutine generates the exchange-correlation potential
  !     in the interstitial region    c.l.fu
  !     including gradient corrections. t.a. 1996.
  !     ******************************************************
CONTAINS
  SUBROUTINE visxcg(ifftd,stars,sym, ifftxc3d, cell,den,&
       xcpot,input, obsolete,noco, vxc,vx,exc)

    !     ******************************************************
    !     instead of visxcor.f: the different exchange-correlation
    !     potentials defined through the key icorr are called through
    !     the driver subroutine vxcallg.f,for the energy density - excallg
    !     subroutines vectorized
    !     ** r.pentcheva 22.01.96
    !     *********************************************************
    !     in case of total = .true. calculates the ex-corr. energy
    !     density
    !     ** r.pentcheva 08.05.96
    !     ******************************************************************

    USE m_grdrsis
    USE m_prpxcfftmap
    USE m_mkgxyz3
    USE m_fft3d
    USE m_fft3dxc
    USE m_types
    IMPLICIT NONE

    CLASS(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_obsolete),INTENT(IN)   :: obsolete
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_potden),INTENT(IN)     :: den
    TYPE(t_potden),INTENT(INOUT)  :: vxc,vx,exc
    
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    !     .. Array Arguments ..

    REAL rhmni ,d_15,sprsv
    !
    !----->  fft  information  for xc potential + energy
    !
    INTEGER, ALLOCATABLE :: igxc_fft(:)
    REAL,    ALLOCATABLE :: gxc_fft(:,:)
    !     ..
    !     .. Local Scalars ..
    INTEGER :: i ,k,js,nt,ifftxc3,idm,jdm,ndm,ig
    COMPLEX :: ci
    REAL    :: rhotot
    INTEGER :: ifftd,ifftxc3d
    TYPE(t_gradients)::grad

    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: fg3(:),cqpw(:,:), ph_wrk(:)
    REAL,    ALLOCATABLE :: bf3(:)
    REAL,    ALLOCATABLE :: rho(:,:),rhd1(:,:,:),rhd2(:,:,:)
    REAL,    ALLOCATABLE :: mx(:),my(:)
    REAL,    ALLOCATABLE :: magmom(:),dmagmom(:,:),ddmagmom(:,:,:) 
    ! 
    REAL, ALLOCATABLE :: v_x(:,:),v_xc(:,:),e_xc(:),vcon(:) 
  
    !     .. unused input (needed for other noco GGA-implementations) ..

    !ta+
    !.....------------------------------------------------------------------
    !------->          abbreviations
    !
    !     ph_wrk: work array containing phase * g_x,gy...... 
    !     den%pw: charge density stored as stars
    !     rho   : charge density stored in real space
    !     v_xc   : exchange-correlation potential in real space
    !     exc   : exchange-correlation energy density in real space
    !     kxc1d  : dimension of the charge density fft box in the pos. domain
    !     kxc2d  : defined in dimens.f program (subroutine apws).1,2,3 indic
    !     kxc3d  ; a_1, a_2, a_3 directions.
    !     kq(i) : i=1,2,3 actual length of the fft-box for which fft is done
    !     nstr  : number of members (arms) of reciprocal lattice (g) vector
    !             of each star.
    !     nxc3_fft: number of stars in the  charge density  fft-box
    !     ng3   : number of 3 dim. stars in the charge density sphere define
    !             by gmax
    !     kmxxc_fft: number of g-vectors forming the nxc3_fft stars in the
    !               charge density or xc-density sphere
    !     kimax : number of g-vectors forming the ng3 stars in the gmax-sphe
    !     ifftxc3d: elements (g-vectors) in the charge density  fft-box
    !     igfft : pointer from the g-sphere (stored as stars) to fft-grid
    !             and     from fft-grid to g-sphere (stored as stars)
    !     pgfft : contains the phases of the g-vectors of sph.
    !     isn   : isn = +1, fft transform for g-space to r-space
    !             isn = -1, vice versa
    !
    !-------------------------------------------------------------------
    !
    !---> set up pointer for backtransformation of from g-vector in
    !     positive domain of xc density fftbox into stars.
    !     also the x,y,z components of the g-vectors are set up to calculate
    !     derivatives.
    !     in principle this can also be done in main program once.
    !     it is done here to save memory.
    !
    ifftd=27*stars%mx1*stars%mx2*stars%mx3
    ifftxc3d = stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft
    ALLOCATE ( igxc_fft(0:ifftxc3d-1),gxc_fft(0:ifftxc3d-1,3) )
    CALL prp_xcfft_map(stars,sym, cell, igxc_fft,gxc_fft)
    !
    ifftxc3=stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft
   
    IF (stars%ng3.GT.stars%ng3) THEN
       WRITE(6,'(/'' stars%ng3.gt.stars%ng3. stars%ng3,stars%ng3='',2i6)') stars%ng3,stars%ng3
       CALL juDFT_error("ng3.gt.n3d",calledby="visxcg")
    ENDIF

    d_15=1.e-15
    !
    ci=CMPLX(0.,1.)
    !
    ! Allocate arrays
    ! ff
    ALLOCATE( bf3(0:ifftd-1),ph_wrk(0:ifftxc3d-1), rho(0:ifftxc3d-1,input%jspins),&
         rhd1(0:ifftxc3d-1,input%jspins,3), rhd2(0:ifftxc3d-1,input%jspins,6) )
    IF (noco%l_noco)  THEN
       ALLOCATE( mx(0:ifftxc3-1),my(0:ifftxc3-1),&
            magmom(0:ifftxc3-1), dmagmom(0:ifftxc3-1,3),ddmagmom(0:ifftxc3-1,3,3) )
    END IF


    !-->     transform charge density to real space

    DO js=1,input%jspins
       CALL fft3dxc(rho(0:,js),bf3, den%pw(:,js), stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
            stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
    END DO

    IF (noco%l_noco) THEN  

       !       for off-diagonal parts the same
       CALL fft3dxc(mx,my, den%pw(:,3), stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
            stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)

       DO i=0,ifftxc3-1 
          rhotot= 0.5*( rho(i,1) + rho(i,2) )
          magmom(i)= SQRT(  (0.5*(rho(i,1)-rho(i,2)))**2 + mx(i)**2 + my(i)**2 )
          rho(i,1)= rhotot+magmom(i)
          rho(i,2)= rhotot-magmom(i)
       END DO

    ENDIF

    IF (.not.xcpot%is_gga()) GOTO 100  

    ! In collinear calculations all derivatives are calculated in g-spce,
    ! in non-collinear calculations the derivatives of |m| are calculated in real space. 

    !-->   for d(rho)/d(x,y,z) = rhd1(:,:,idm) (idm=1,2,3).
    !
    !         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z).

    ALLOCATE(cqpw(stars%ng3,input%jspins))

    DO js= 1,input%jspins
       DO i = 1,stars%ng3
          cqpw(i,js)= ci*den%pw(i,js)
       END DO
    END DO

    DO idm=1,3

       DO ig = 0 , stars%kmxxc_fft - 1
          ph_wrk(ig) = stars%pgfft(ig) * gxc_fft(ig,idm)
       END DO

       DO js=1,input%jspins
          CALL fft3dxc(rhd1(0:,js,idm),bf3, cqpw(:,js), stars%kxc1_fft,stars%kxc2_fft,&
               stars%kxc3_fft,stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,ph_wrk,stars%nstr)
       END DO

    END DO

    IF (noco%l_noco) THEN

       CALL grdrsis(magmom,cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,obsolete, dmagmom )

       DO i=0,ifftxc3-1
          DO idm=1,3
             rhotot= rhd1(i,1,idm)/2.+rhd1(i,2,idm)/2.
             rhd1(i,1,idm)= rhotot+dmagmom(i,idm) 
             rhd1(i,2,idm)= rhotot-dmagmom(i,idm) 
          END DO
       END DO

    END IF

   
    !-->   for dd(rho)/d(xx,xy,yy,zx,yz,zz) = rhd2(:,:,idm) (idm=1,2,3,4,5,6)
    !
    !         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z) * g_(x,y,z)

    DO i = 1,stars%ng3
       DO js=1,input%jspins 
          cqpw(i,js)= -den%pw(i,js)
       END DO
    END DO

    ndm = 0
    DO idm = 1,3
       DO jdm = 1,idm
          ndm = ndm + 1

          DO ig = 0 , stars%kmxxc_fft-1
             ph_wrk(ig) = stars%pgfft(ig)*gxc_fft(ig,idm)*gxc_fft(ig,jdm)
          ENDDO

          DO js=1,input%jspins
             CALL fft3dxc(rhd2(0:,js,ndm),bf3, cqpw(:,js), stars%kxc1_fft,stars%kxc2_fft,&
                  stars%kxc3_fft,stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,ph_wrk,stars%nstr)
          END DO
       END DO ! jdm 
    END DO   ! idm 

    DEALLOCATE(cqpw)

    IF (noco%l_noco) THEN

       DO idm = 1,3
          CALL grdrsis(dmagmom(0,idm),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,obsolete, ddmagmom(0,1,idm) )
       END DO

       ndm= 0
       DO idm = 1,3
          DO jdm = 1,idm
             ndm = ndm + 1  

             DO i=0,ifftxc3-1
                rhotot= rhd2(i,1,ndm)/2.+rhd2(i,2,ndm)/2.
                rhd2(i,1,ndm)= rhotot + ( ddmagmom(i,jdm,idm) + ddmagmom(i,idm,jdm) )/2. 
                rhd2(i,2,ndm)= rhotot - ( ddmagmom(i,jdm,idm) + ddmagmom(i,idm,jdm) )/2. 
             END DO

          ENDDO !jdm
       ENDDO   !idm 

    END IF

100 CONTINUE


    DEALLOCATE ( ph_wrk )
    IF (noco%l_noco) THEN 
       DEALLOCATE(mx,my,magmom,dmagmom,ddmagmom) 
    END IF
    !
    DO js=1,input%jspins 
       DO i=0,ifftxc3-1
          rho(i,js)=MAX(rho(i,js),d_15)
       ENDDO
    END DO
    bf3=0.0
    ! allocate the other arrays 
    !
    CALL xcpot%alloc_gradients(ifftxc3d,input%jspins,grad)
 
    !
    !     calculate the quantities such as abs(grad(rho)),.. used in
    !     evaluating the gradient contributions to potential and energy.
    !
    CALL mkgxyz3 (ifftxc3d,input%jspins,ifftxc3,input%jspins,rho, rhd1(0,1,1),rhd1(0,1,2),rhd1(0,1,3),&
         rhd2(0,1,1),rhd2(0,1,3),rhd2(0,1,6), rhd2(0,1,5),rhd2(0,1,4),rhd2(0,1,2), grad)

    DEALLOCATE ( rhd1,rhd2 )
    ALLOCATE ( v_xc(0:ifftxc3d-1,input%jspins) )
    ALLOCATE ( v_x (0:ifftxc3d-1,input%jspins) )
    !
    !     calculate the exchange-correlation potential in  real space
    !
    nt=ifftxc3
    !
    !      rhmni: rho_minimum_interstitial.

    rhmni=10.e+10

    DO js=1,input%jspins
       DO i=0,ifftxc3-1
          rho(i,js)=MAX(rho(i,js),d_15)
          rhmni=MIN(rhmni,rho(i,js))
       ENDDO
    ENDDO

  
    IF (rhmni.LT.obsolete%chng) THEN
       WRITE(6,'(/'' rhmn.lt.obsolete%chng in visxc. rhmn,obsolete%chng='',2d9.2)') rhmni,obsolete%chng
       !          CALL juDFT_error("visxcg: rhmn.lt.chng",calledby="visxcg")
    ENDIF
    CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)
    !
    !----> back fft to g space
    !----> perform back  fft transform: v_xc(r) --> vxc(star)
    !
    ALLOCATE(fg3(stars%ng3))

    DO js = 1,input%jspins
       bf3=0.0
       CALL fft3dxc(v_xc(0:,js),bf3, fg3, stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
            stars%nxc3_fft,stars%kmxxc_fft,-1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
       !
       DO k = 1,stars%nxc3_fft
          vxc%pw(k,js) = vxc%pw(k,js) + fg3(k)
       ENDDO

       !
       !====>  INCLUDING TOTAL ENERGY
       !
       IF (input%total) THEN
          !
          !----> Perform fft transform: v_xc(star) --> vxc(r) 
          !     !Use large fft mesh for convolution
          !
          fg3(stars%nxc3_fft+1:)=0.0
          ALLOCATE ( vcon(0:ifftd-1) )
          CALL fft3d(vcon(0),bf3, fg3, stars,+1)
          !
          !----> Convolute with step function
          !
          DO i=0,ifftd-1
             vcon(i)=stars%ufft(i)*vcon(i)
          ENDDO
          bf3=0.0
          CALL fft3d(vcon(0),bf3, fg3, stars,-1,.FALSE.)
          DEALLOCATE ( vcon )
          !
          !----> add to warped coulomb potential
          !
          DO k = 1,stars%ng3
             vxc%pw_w(k,js) = vxc%pw_w(k,js) + fg3(k)
          ENDDO

       ENDIF

       bf3=0.0
       CALL fft3dxc(v_x(0:,js),bf3, fg3, stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft, &
            stars%nxc3_fft,stars%kmxxc_fft,-1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
       !
       DO k = 1,stars%nxc3_fft
          vx%pw(k,js) = vx%pw(k,js) + fg3(k)
       ENDDO

       !
       !====>   INCLUDING TOTAL ENERGY
       !
       IF (input%total) THEN
          !
          !---->  Perform fft transform: v_xc(star) --> vxc(r) 
          !       !Use large fft mesh for convolution
          !

          fg3(stars%nxc3_fft+1:)=0.0
          ALLOCATE ( vcon(0:ifftd-1) )
          CALL fft3d(vcon(0),bf3, fg3, stars,+1)
          !
          !----> Convolute with step function
          !
          DO i=0,ifftd-1
             vcon(i)=stars%ufft(i)*vcon(i)
          ENDDO
          bf3=0.0
          CALL fft3d(vcon(0),bf3, fg3, stars,-1,.FALSE.)
          DEALLOCATE ( vcon )
          !
          !----> add to warped exchange-potential
          !
          DO k = 1,stars%ng3
             vx%pw_w(k,js) = vx%pw_w(k,js) + fg3(k)
          ENDDO

       ENDIF

    ENDDO
    DEALLOCATE ( v_x,v_xc )
    !
    !     calculate the ex.-cor energy density in real space
    !
    IF (ALLOCATED(exc%pw_w)) THEN
       ALLOCATE ( e_xc(0:ifftxc3d-1) )
       CALL xcpot%get_exc(input%jspins,rho,e_xc,grad)
       !
       !---->   perform back  fft transform: exc(r) --> exc(star)
       !
       bf3=0.0
       CALL fft3dxc(e_xc,bf3, fg3, stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
            stars%nxc3_fft,stars%kmxxc_fft,-1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
       DEALLOCATE ( e_xc )
       !
       !---->   Perform fft transform: exc(star) --> exc(r) 
       !        !Use large fft mesh for convolution
       !
       fg3(stars%nxc3_fft+1:)=0.0
       bf3=0.0
       ALLOCATE ( vcon(0:ifftd-1) )
       CALL fft3d(vcon,bf3,fg3, stars,+1)

       DO i=0,ifftd-1
          vcon(i)=stars%ufft(i)*vcon(i)
       ENDDO
       !
       !         ---> back fft to g space
       !
       bf3=0.0
       CALL fft3d(vcon,bf3,exc%pw_w(:,1), stars,-1,.FALSE.)
       DEALLOCATE ( vcon )
       !
    ENDIF

    DEALLOCATE(fg3)
    DEALLOCATE ( bf3,rho,igxc_fft,gxc_fft )
  



  END SUBROUTINE visxcg
END MODULE m_visxcg
