MODULE m_vvacxc
  !     ********************************************************************
  !     calculates 2-dim star function coefficients of exchange-correlation*
  !     potential in the vacuum regions  and adds them to the corresponding*
  !     coeffs of the coulomb potential            c.l.fu, r.podloucky     *
  !     ********************************************************************
CONTAINS
  SUBROUTINE vvacxc(ifftd2,stars,vacuum,xcpot,input,noco,den,vxc,exc)

    !     ********************************************************************
    !     instead of vvacxcor.f: the different exchange-correlation 
    !     potentials defined through the key icorr are called through 
    !     the driver subroutine vxcall.f, subroutines vectorized
    !     in case of TOTAL = .TRUE. calculates the ex.-corr. energy
    !     density through the driver subroutine excall.f
    !     ** r.pentcheva 08.05.96
    !     ********************************************************************

    USE m_xcall, ONLY : vxcall,excall
    USE m_fft2d
    USE m_types
    IMPLICIT NONE
    TYPE(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_potden),INTENT(IN)    :: den
    TYPE(t_potden),INTENT(INOUT) :: vxc,exc
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ifftd2
    !     ..
 
    !     ..
    !     .. Local Scalars ..
    INTEGER :: k,js,nt,irec2,nmzdiff,ivac,ip,i 
    REAL    :: rhti
    REAL    :: chdens,magmom 
    !     ..
    !     .. Local Arrays ..
    COMPLEX :: fgxy(stars%ng2-1)
    REAL    :: af2(0:ifftd2-1,input%jspins),bf2(0:ifftd2-1),fgz
    REAL,ALLOCATABLE :: mx(:),my(:) 
    !     warping region
    REAL    :: v_xc(0:ifftd2-1,input%jspins),e_xc(0:ifftd2-1)
    REAL    :: v_x (0:ifftd2-1,input%jspins)
    !     beyond warping region
    REAL    :: vxcz(vacuum%nmzd,input%jspins)
    REAL    :: vxz (vacuum%nmzd,input%jspins)

    IF (noco%l_noco) THEN
       ALLOCATE (mx(0:ifftd2-1),my(0:ifftd2-1)) 
    END IF

    nt = ifftd2
    rhti = 0.
    !
    !     the charge density in vacuum is expanded in 2-dim stars on a mesh 
    !     in z-direction . The G||.ne.zero-components expand from 1 to nmzxy
    !     the G||.eq.zero-components expand from 1 to nmz
    !     first we calculate vxc in the warping region
    !
    DO  ivac = 1,vacuum%nvac
       DO ip = 1,vacuum%nmzxy
          !
          !         transform charge density to real space: 2-dim FFT
          !
          DO js = 1,input%jspins
             CALL fft2d(&
                  &               stars,&
                  &               af2(0,js),bf2,&
                  &               den%vacz(ip,ivac,js),rhti,den%vacxy(ip,1,ivac,js),&
                  &               vacuum%nmzxyd,+1)
          END DO

          IF (noco%l_noco) THEN 

             CALL fft2d(&
                  &               stars,&
                  &               mx,my, &
                  &               den%vacz(ip,ivac,3),den%vacz(ip,ivac,4),&
                  &               den%vacxy(ip,1,ivac,3),&
                  &               vacuum%nmzxyd,1)
             DO i=0,9*stars%mx1*stars%mx2-1 
                chdens= (af2(i,1)+af2(i,2))/2.  
                magmom= mx(i)**2 + my(i)**2 +&
                     &                ((af2(i,1)-af2(i,2))/2.)**2 
                magmom= SQRT(magmom) 
                af2(i,1)= chdens + magmom 
                af2(i,2)= chdens - magmom
             END DO

          END IF
          !
          !         calculate the exchange-correlation potential in  real space
          !
          CALL vxcall&
               &               (6,xcpot,input%jspins,&
               &                ifftd2,nt,af2,&
               &                v_x,v_xc)   

          DO  js = 1,input%jspins
             !
             !            ----> 2-d back fft to g space
             !
             bf2=0.0
             CALL fft2d(&
                  &                 stars,&
                  &                 v_xc(0,js),bf2,&
                  &                 fgz,rhti,fgxy,&
                  &                 1,-1)
             !
             !            ----> and add vxc to coulomb potential
             !            the G||.eq.zero component is added to vz
             !
             vxc%vacz(ip,ivac,js) = fgz + vxc%vacz(ip,ivac,js)
             !
             !            the G||.ne.zero components are added to vxc%vacxy
             !
             DO irec2 = 1,stars%ng2-1
                vxc%vacxy(ip,irec2,ivac,js) = vxc%vacxy(ip,irec2,ivac,js) +&
                     &                                              fgxy(irec2)
             ENDDO
          ENDDO
          !
          !i        calculate the exchange-correlation energy density in  real space
          !
          IF (input%total) THEN   
             CALL excall&
                  &                 (6,xcpot,input%jspins,&
                  &                  ifftd2,nt,af2,&
                  &                  e_xc)   
             !
             !     ----> 2-d back fft to g space
             !
             bf2=0.0
             CALL fft2d(&
                  &                 stars,&
                  &                 e_xc,bf2,&
                  &                 exc%vacz(ip,ivac,1),rhti,exc%vacxy(ip,1,ivac,1),&
                  &                 vacuum%nmzxyd,-1)
          ENDIF

       ENDDO
       !
       !        calculate vxc for z now beyond warping region 
       !
       nmzdiff = vacuum%nmz - vacuum%nmzxy

       DO k=1,nmzdiff

          DO js=1,input%jspins
             af2(k-1,js) = den%vacz(vacuum%nmzxy+k,ivac,js)
          ENDDO

          IF (noco%l_noco) THEN

             mx(0)= den%vacz(vacuum%nmzxy+k,ivac,3)
             my(0)= den%vacz(vacuum%nmzxy+k,ivac,4)
             chdens= (af2(k-1,1)+af2(k-1,2))/2.
             magmom= mx(0)**2 + my(0)**2 +&
                  &               ((af2(k-1,1)-af2(k-1,2))/2.)**2
             magmom= SQRT(magmom)
             af2(k-1,1)= chdens + magmom
             af2(k-1,2)= chdens - magmom

          END IF

       ENDDO

       CALL vxcall&
            &              (6,xcpot,input%jspins,&
            &               vacuum%nmzd,nmzdiff,af2,&
            &               vxz,vxcz)
       !+gu
       DO  js=1,input%jspins
          DO k=vacuum%nmzxy+1,vacuum%nmz
             vxc%vacz(k,ivac,js) = vxc%vacz(k,ivac,js) + vxcz(k-vacuum%nmzxy,js)
          ENDDO
       ENDDO
       !
       WRITE (6,FMT=8020) ivac, (vxc%vacz(vacuum%nmz,ivac,js),js=1,input%jspins)
       WRITE (16,FMT=8020) ivac, (vxc%vacz(vacuum%nmz,ivac,js),js=1,input%jspins)
8020   FORMAT (/,5x,'vacuum zero for vacuum',i3,' = ',2f10.5)
       !
       !        calculate the ex.-corr. energy density now beyond warping region
       !
       IF (input%total) THEN
          CALL excall&
               &                   (6,xcpot,input%jspins,&
               &                    vacuum%nmzd,nmzdiff,af2,&
               &                    exc%vacz(vacuum%nmzxy+1,ivac,1))
       END IF
    ENDDO
    IF (noco%l_noco) THEN 
       DEALLOCATE (mx,my)
    END IF

  END SUBROUTINE vvacxc
END MODULE m_vvacxc
