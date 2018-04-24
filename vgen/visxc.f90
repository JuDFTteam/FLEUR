      MODULE m_visxc
!     ******************************************************
!     subroutine generates the exchange-correlation potential
!     in the interstitial region    c.l.fu
!     ******************************************************
      CONTAINS
      SUBROUTINE visxc(ifftd,stars,noco,xcpot,input,den,vxc,vx,exc)

!     ******************************************************
!     instead of visxcor.f: the different exchange-correlation 
!     potentials defined through the key icorr are called through 
!     the driver subroutine vxcall.f,for the energy density - excall
!     subroutines vectorized
!     in case of TOTAL = .TRUE. calculates the ex.-corr. energy density
!     ** r.pentcheva 08.05.96
!     ********************************************************************
      USE m_types
      USE m_fft3d
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)         :: ifftd
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_noco),INTENT(IN)      :: noco
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_potden),INTENT(IN)    :: den
      TYPE(t_potden),INTENT(INOUT) :: vxc,exc,vx
      
           
!     ..
!     .. Local Scalars ..
      INTEGER i,k,js,nt
      REAL    chdens,magmom
!     ..
!     .. Local Arrays
      COMPLEX fg3(stars%ng3)
      REAL, ALLOCATABLE :: mx(:),my(:)
      REAL, ALLOCATABLE :: e_xc(:),vcon(:),v_xc(:,:),v_x(:,:)
      REAL, ALLOCATABLE :: af3(:,:),bf3(:)
!
!     ---> allocate arrays
!
      ALLOCATE ( e_xc(0:ifftd-1),vcon(0:ifftd-1),v_xc(0:ifftd-1,input%jspins),&
           af3(0:ifftd-1,input%jspins),bf3(0:ifftd-1) )

      ALLOCATE( v_x(0:ifftd-1,input%jspins) )
!
!     transform charge density to real space
!
      DO js = 1,input%jspins
         CALL fft3d(af3(0,js),bf3, den%pw(1,js), stars,+1)
      ENDDO

      IF (noco%l_noco) THEN 

        ALLOCATE (mx(0:ifftd-1),my(0:ifftd-1))

        CALL fft3d(mx,my, den%pw(:,3), stars,+1)
        DO i=0,27*stars%mx1*stars%mx2*stars%mx3-1
          chdens= (af3(i,1)+af3(i,2))/2.
          magmom= mx(i)**2 + my(i)**2 + ((af3(i,1)-af3(i,2))/2.)**2
          magmom= SQRT(magmom)
          af3(i,1)= chdens + magmom
          af3(i,2)= chdens - magmom
        END DO

        DEALLOCATE (mx,my)

      END IF

!
!     calculate the exchange-correlation potential in  real space
!
      nt=ifftd
      CALL xcpot%get_vxc(input%jspins,af3,v_xc,v_x)
      
!
!     ---> back fft to g space and add to coulomb potential for file nrp
!

      IF (ALLOCATED(exc%pw_w)) THEN

         DO js = 1,input%jspins

           DO i=0,ifftd-1
             vcon(i)=stars%ufft(i)*v_xc(i,js)
             bf3(i)=0.0
           ENDDO
           CALL fft3d(vcon,bf3, fg3, stars,-1)
           fg3=fg3*stars%nstr

           DO k = 1,stars%ng3
              vxc%pw_w(k,js) = vxc%pw_w(k,js) + fg3(k)
           ENDDO   

           DO i=0,ifftd-1
             vcon(i)=stars%ufft(i)*v_x(i,js)
             bf3(i)=0.0
           ENDDO
           CALL fft3d(vcon,bf3, fg3, stars,-1)
           fg3=fg3*stars%nstr

           DO k = 1,stars%ng3
              vx%pw_w(k,js) = vx%pw_w(k,js) + fg3(k)
           ENDDO   

         ENDDO
!
!     calculate the ex.-cor energy density in real space
         !
         CALL xcpot%get_exc(input%jspins,af3,e_xc)

         DO i=0,ifftd-1
           vcon(i)=stars%ufft(i)*e_xc(i)
           bf3(i)=0.0
         ENDDO
!
!         ---> back fft to g space
!
         CALL fft3d(vcon,bf3, exc%pw_w(:,1), stars,-1)
         exc%pw_w(:,1)=exc%pw_w(:,1)*stars%nstr
!
      END IF ! input%total

      DO js = 1,input%jspins
         bf3(0:ifftd-1)=0.0
         CALL fft3d(v_xc(0,js),bf3, fg3, stars,-1)
         DO k = 1,stars%ng3
            vxc%pw(k,js) = vxc%pw(k,js) + fg3(k)
         ENDDO
         
         CALL fft3d(v_x(0,js),bf3, fg3, stars,-1)
         DO k = 1,stars%ng3
            vx%pw(k,js) = vx%pw(k,js) + fg3(k)
         ENDDO      
      
      ENDDO

      END SUBROUTINE visxc
      END MODULE m_visxc
