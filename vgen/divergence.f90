!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_divergence
   USE m_types
   USE m_juDFT
   PRIVATE
   PUBLIC :: divergence, vac_grad, divpotgrad

CONTAINS
   SUBROUTINE divergence(input,stars,atoms,sphhar,vacuum,sym,cell,noco,bxc,div)
      USE m_lattHarmsSphHarmsConv
      USE m_gradYlm
      USE m_constants

      !--------------------------------------------------------------------------
      ! Use the interstitial/vacuum divergence subroutine and an external MT-gra-
      ! dient routine from juPhon to assemble the divergence of a field into a
      ! t_potden variable. The MT-gradient is first calculated in sperical har-
      ! monics coefficients.
      !--------------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(t_input),                INTENT(IN)    :: input
      TYPE(t_stars),                INTENT(IN)    :: stars
      TYPE(t_atoms),                INTENT(IN)    :: atoms
      TYPE(t_sphhar),               INTENT(IN)    :: sphhar
      TYPE(t_vacuum),               INTENT(IN)    :: vacuum
      TYPE(t_sym),                  INTENT(IN)    :: sym
      TYPE(t_cell),                 INTENT(IN)    :: cell
      TYPE(t_noco),                 INTENT(IN)    :: noco
      TYPE(t_potden), DIMENSION(3), INTENT(INOUT) :: bxc
      TYPE(t_potden),               INTENT(INOUT) :: div

      TYPE(t_potden), DIMENSION(3)                :: grad

      INTEGER :: i,iType,indmax, lh
      COMPLEX, ALLOCATABLE :: flm(:,:,:),grsflm1(:,:,:,:),grsflm2(:,:,:,:),grsflm3(:,:,:,:),divflm(:,:,:) ! (iR,lm,n[,x,i])

      CALL timestart("MT divergence")
      indmax=(atoms%lmaxd+1)**2

      ALLOCATE(flm(atoms%jmtd,indmax,atoms%ntype))
      ALLOCATE(divflm(atoms%jmtd,indmax,atoms%ntype))

      DO i=1,3
         DO iType=1, atoms%ntype
            CALL lattHarmsRepToSphHarms(sym, atoms, sphhar, iType, bxc(i)%mt(:,:,iType,1), flm(:,:,iType))
         END DO
         IF (i==1) THEN
            CALL gradYlm(atoms,flm,grsflm1)
         ELSE IF (i==2) THEN
            CALL gradYlm(atoms,flm,grsflm2)
         ELSE
            CALL gradYlm(atoms,flm,grsflm3)
         END IF
      END DO

      DEALLOCATE(flm)

      CALL divYlm(grsflm1(:,:indmax,:,:),grsflm2(:,:indmax,:,:),grsflm3(:,:indmax,:,:), divflm)

      DO iType=1, atoms%ntype
         CALL sphHarmsRepToLattHarms(sym, atoms, sphhar, iType, divflm(:,1:indmax,iType), div%mt(:,0:,iType,1))
      END DO

      DEALLOCATE(divflm,grsflm1,grsflm2,grsflm3)

      CALL timestop("MT divergence")

      CALL timestart("PW divergence")

      div%pw(:,1)=CMPLX(0.0,0.0)

      DO i=1,3
         div%pw(:,1)=div%pw(:,1)+ImagUnit*(cell%bmat(i,1)*stars%kv3(1,:)+cell%bmat(i,2)*stars%kv3(2,:)+cell%bmat(i,3)*stars%kv3(3,:))*bxc(i)%pw(:,1)
      END DO

      CALL timestop("PW divergence")

      IF (input%film) THEN
         CALL timestart("Vac divergence")
         div%vacxy=CMPLX(0.0,0.0)
         div%vacz=0.0
         CALL vac_grad(vacuum,stars,bxc(1),grad,9*stars%mx1*stars%mx2)
         div%vacxy=div%vacxy+grad(1)%vacxy
         div%vacz=div%vacz+grad(1)%vacz
         CALL vac_grad(vacuum,stars,bxc(2),grad,9*stars%mx1*stars%mx2)
         div%vacxy=div%vacxy+grad(2)%vacxy
         div%vacz=div%vacz+grad(2)%vacz
         CALL vac_grad(vacuum,stars,bxc(3),grad,9*stars%mx1*stars%mx2)
         div%vacxy=div%vacxy+grad(3)%vacxy
         div%vacz=div%vacz+grad(3)%vacz
         CALL timestop("Vac divergence")
      END IF


   END SUBROUTINE divergence

   SUBROUTINE vac_grad(vacuum,stars,den,grad,ifftd2)

      USE m_constants
      USE m_grdchlh
      USE m_fft2d
      USE m_types

      IMPLICIT NONE
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_potden),INTENT(IN)    :: den
      TYPE(t_potden),INTENT(INOUT),DIMENSION(3) :: grad
      !     ..
      !     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ifftd2

      !     ..
      !     .. Local Scalars ..
      INTEGER :: js,nt,i,iq,irec2,nmz0,nmzdiff,ivac,ip
      REAL    :: rhti,zro,d_15
      !     ..
      !     .. Local Arrays ..
      REAL :: fgz(3)
      REAL, ALLOCATABLE :: af2(:),bf2(:)
      REAL, ALLOCATABLE :: rhdx(:),rhdy(:),rhdz(:)
      REAL, ALLOCATABLE :: rhtdz(:),rhtdzz(:)
      REAL, ALLOCATABLE :: rxydzr(:),rxydzi(:)
      REAL, ALLOCATABLE :: rxydzzr(:),rxydzzi(:),rhtxyr(:),rhtxyi(:)
      REAL, ALLOCATABLE :: rhtxc(:,:),dummy(:)
      COMPLEX, ALLOCATABLE :: fgxy(:,:),rxydz(:,:),rxydzz(:),cqpw(:)

      d_15     = 1.e-15
      zro      = 0.0
      nt       = ifftd2

      ALLOCATE (rxydz(vacuum%nmzxy,stars%ng2-1))
      ALLOCATE (rhtdz(vacuum%nmzd),rhtdzz(vacuum%nmzd))

      DO ivac=1,vacuum%nvac

         ! the charge density in vacuum is expanded in 2-dim stars on a mesh
         ! in z-direction. the g||.ne.zero-components expand from 1 to nmzxy
         ! the g||.eq.zero-components expand from 1 to nmz
         ! first we calculate vxc in the warping region

         !
         ! calculate first (rhtdz) & second (rhtdzz) derivative of den%vacz(1:nmz)
         !

         CALL grdchlh(vacuum%delz,den%vacz(1:vacuum%nmz,ivac,1),&
                     rhtdz,rhtdzz)
         ALLOCATE ( rhtxyr(vacuum%nmzxy), rhtxyi(vacuum%nmzxy),dummy(vacuum%nmzxy) )
         ALLOCATE ( rxydzr(vacuum%nmzxy), rxydzi(vacuum%nmzxy) )
         ALLOCATE ( rxydzzr(vacuum%nmzxy),rxydzzi(vacuum%nmzxy) )

         DO iq = 1, stars%ng2-1
         !
         ! calculate first (rxydz) & second (rxydzz) derivative of den%vacxy:
         !
            DO ip=1,vacuum%nmzxy
               rhtxyr(ip)=den%vacxy(ip,iq,ivac,1)
            ENDDO
            CALL grdchlh(vacuum%delz,rhtxyr(:vacuum%nmzxy), rxydzr,rxydzzr)

            DO ip=1,vacuum%nmzxy
               rhtxyi(ip)=aimag(den%vacxy(ip,iq,ivac,js))
            ENDDO

            CALL grdchlh(vacuum%delz,rhtxyi(:vacuum%nmzxy), rxydzi,rxydzzi)

            DO ip=1,vacuum%nmzxy
               rxydz(ip,iq)=cmplx(rxydzr(ip),rxydzi(ip))
            ENDDO

         ENDDO ! loop over 2D stars (iq)

         DEALLOCATE ( rhtxyr,rhtxyi,rxydzr,rxydzi,rxydzzr,rxydzzi )
         DEALLOCATE ( dummy )

         ALLOCATE ( rhdx(0:ifftd2-1),rhdy(0:ifftd2-1) )
         ALLOCATE ( rhdz(0:ifftd2-1))

         ALLOCATE ( cqpw(stars%ng2-1),af2(0:ifftd2-1) )
         ALLOCATE ( fgxy(stars%ng2-1,3),bf2(0:ifftd2-1) )

         af2=0.0
         DO ip = 1,vacuum%nmzxy
            ! loop over warping region

            ! Transform charge and magnetization to real-space.

            CALL fft2d(stars, af2(0),bf2, den%vacz(ip,ivac,1),0.,&
                       den%vacxy(ip,1,ivac,1), vacuum%nmzxyd,+1)

            ! calculate derivatives with respect to x,y in g-space
            ! and transform them to real-space.

            DO iq=1,stars%ng2-1
               cqpw(iq)=ImagUnit*den%vacxy(ip,iq,ivac,js)
            ENDDO

            rhti = 0.0
            ! d(rho)/atoms%dx is obtained by a FFT of i*gx*den%vacxy
            ! (den%vacz is set to zero and gx is included in
            ! dn/atoms =  FFT(0,i*gx*den%vacxy)
            CALL fft2d(stars, rhdx(0),bf2, zro,rhti,cqpw, 1,+1,stars%ft2_gfx)

            rhti = 0.0
            CALL fft2d(    &               ! dn/dy =  FFT(0,i*gy*den%vacxy)&
                        stars, rhdy(0),bf2, zro,rhti,cqpw, 1,+1,stars%ft2_gfy)

            rhti = 0.0
            CALL fft2d(     &              ! dn/dz = FFT(rhtdz,rxydz)&
                      stars, rhdz(0),bf2, rhtdz(ip),rhti,rxydz(ip,1), vacuum%nmzxyd,+1)

            !
            ! set minimal value of af2 to 1.0e-15
            !

            !af2=max(af2,10e-13)

            where (af2<d_15) af2=d_15

            !
            !           ----> 2-d back fft to g space
            !
            bf2=0.0
            CALL fft2d(stars, rhdx,bf2, fgz(1),rhti,fgxy(:,1), 1,-1)
            CALL fft2d(stars, rhdy,bf2, fgz(2),rhti,fgxy(:,2), 1,-1)
            CALL fft2d(stars, rhdz,bf2, fgz(3),rhti,fgxy(:,3), 1,-1)

            ! the g||.eq.zero component is added to grad%vacz
            !
            grad(1)%vacz(ip,ivac,1) = fgz(1) + grad(1)%vacz(ip,ivac,1)
            grad(2)%vacz(ip,ivac,1) = fgz(2) + grad(2)%vacz(ip,ivac,1)
            grad(3)%vacz(ip,ivac,1) = fgz(3) + grad(3)%vacz(ip,ivac,1)
            !
            ! the g||.ne.zero components are added to grad%vacxy
            !
            DO irec2 = 1,stars%ng2-1
               grad(1)%vacxy(ip,irec2,ivac,1)=grad(1)%vacxy(ip,irec2,ivac,1)+fgxy(irec2,1)
               grad(2)%vacxy(ip,irec2,ivac,1)=grad(2)%vacxy(ip,irec2,ivac,1)+fgxy(irec2,2)
               grad(3)%vacxy(ip,irec2,ivac,1)=grad(3)%vacxy(ip,irec2,ivac,1)+fgxy(irec2,3)
            ENDDO

         END DO ! ip=1,vacuum%nmzxy
         DEALLOCATE ( rhdx,rhdy,rhdz)
         DEALLOCATE ( cqpw,fgxy)

         ! now treat the non-warping region

         nmzdiff = vacuum%nmz - vacuum%nmzxy

         ! The non-warping region runs from nmzxy+1 to nmz.
         ! The values from nmz0 to nmzxy are taken into account in order
         ! to get the real-space derivative smooth around nmzxy+1.

         nmz0= vacuum%nmzxy+1+(6/2)-6
         IF (nmz0 <= 0) THEN ! usually vacuum%nmzxy>6
            nmz0= 1
         END IF

         DEALLOCATE ( af2)

         DO ip = vacuum%nmzxy + 1,vacuum%nmz
            grad(3)%vacz(ip,ivac,1) = grad(3)%vacz(ip,ivac,1) + rhtdz(ip-vacuum%nmzxy)
         ENDDO

         DEALLOCATE ( bf2)

      ENDDO    ! loop over vacua (ivac)

   END SUBROUTINE vac_grad

   SUBROUTINE divpotgrad(input,stars,atoms,sphhar,vacuum,sym,cell,noco,pot,grad)

      USE m_types
      USE m_lattHarmsSphHarmsConv
      USE m_gradYlm
      USE m_constants

      !--------------------------------------------------------------------------
      ! Use the interstitial/vacuum gradient subroutine and an external MT-gra-
      ! dient routine from juPhon to assemble the gradient of a potenital into a
      ! t_potden variable. The MT-gradient is first calculated in sperical har-
      ! monics coefficients.
      !--------------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(t_input), INTENT(IN)                   :: input
      TYPE(t_stars),INTENT(IN)                    :: stars
      TYPE(t_atoms), INTENT(IN)                   :: atoms
      TYPE(t_sphhar), INTENT(IN)                  :: sphhar
      TYPE(t_vacuum),INTENT(IN)                   :: vacuum
      TYPE(t_sym), INTENT(IN)                     :: sym
      TYPE(t_cell),INTENT(IN)                     :: cell
      TYPE(t_noco), INTENT(IN)                    :: noco
      TYPE(t_potden), INTENT(IN)                  :: pot
      TYPE(t_potden), dimension(3), INTENT(INOUT) :: grad

      TYPE(t_potden)                              :: denloc
      INTEGER :: i,iType,indmax,lh,lhmax
      COMPLEX, ALLOCATABLE :: flm(:,:,:),grsflm(:,:,:,:) ! (iR,lm,n[,x,i])

      CALL timestart("MT potential gradient")
      indmax=(atoms%lmaxd+1)**2

      ALLOCATE(flm(atoms%jmtd,indmax,atoms%ntype))

      denloc=pot

      DO iType=1,atoms%ntype
         lhmax=sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:iType - 1)) + 1))
         DO lh=0, lhmax
            denloc%mt(:,lh,iType,1) = denloc%mt(:,lh,iType,1)*atoms%rmsh(:, iType)**2
         END DO ! lh
         CALL lattHarmsRepToSphHarms(sym, atoms, sphhar, iType, denloc%mt(:,:,iType,1), flm(:,:,iType))
      END DO

      CALL gradYlm(atoms,flm,grsflm)

      DEALLOCATE(flm)

      DO i=1,3
         DO iType=1,atoms%ntype
            CALL sphHarmsRepToLattHarms(sym, atoms, sphhar, iType, grsflm(:,1:indmax,iType,i)/(4.0*pi_const), grad(i)%mt(:,0:,iType,1))
         END DO
      END DO

      DEALLOCATE(grsflm)

      CALL timestop("MT potential gradient")

      CALL timestart("PW potential gradient")

      DO i=1,3
         grad(i)%pw(:,1)=ImagUnit*(cell%bmat(i,1)*stars%kv3(1,:)+cell%bmat(i,2)*stars%kv3(2,:)+cell%bmat(i,3)*stars%kv3(3,:))*pot%pw(:,1)/(4.0*pi_const)
      END DO


      CALL timestop("PW potential gradient")

      IF (input%film) THEN
         CALL timestart("Vac potential gradient")
         CALL vac_grad(vacuum,stars,pot,grad,9*stars%mx1*stars%mx2)
         CALL timestart("Vac potential gradient")
      END IF

   END SUBROUTINE divpotgrad
END MODULE m_divergence
