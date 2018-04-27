MODULE m_qintsl
  USE m_juDFT
CONTAINS
  SUBROUTINE q_int_sl(isp,ikpt,stars,atoms,sym,cell,ne,lapw,slab,oneD,zMat)          
    !     *******************************************************
    !     calculate the charge of the En(k) state 
    !     in the interstitial region of each leyer
    !                                             Yu.M. Koroteev
    !             From pwden_old.F and pwint.F by  c.l.fu
    !     *******************************************************
#include"cpp_double.h"
    USE m_pwintsl
    USE m_types
    IMPLICIT NONE

    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_stars),INTENT(IN)  :: stars
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_zMat),INTENT(IN)   :: zMat
    TYPE(t_slab),INTENT(INOUT):: slab
    !
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne,isp,ikpt
    !     ..
    !     .. Local Scalars ..
    REAL q1,zsl1,zsl2,qi,volsli,volintsli
    INTEGER i ,indp,ix1,iy1,iz1,j,n,ns,ind
    COMPLEX x,phase,phasep
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: stfunint(:,:),z_z(:)

    !     ..
    IF (oneD%odi%d1) CALL juDFT_error("well, does not work with 1D. Not clear how to define a layer.",calledby ="q_int_sl")
    !
    !     calculate the star function expansion coefficients of
    !     the plane wave charge density for each En(k)
    !    
    !     ----> g=0 star
    !
    ALLOCATE ( stfunint(stars%ng3,slab%nsl), z_z(stars%ng3) ) 
    !
    !  -----> calculate the integrals of star functions over
    !                     the layer interstitial
    !
    DO i = 1,slab%nsl
       zsl1 = slab%zsl(1,i)
       zsl2 = slab%zsl(2,i)
       volsli = slab%volsl(i)
       volintsli = slab%volintsl(i)
       DO j = 1,stars%ng3
          CALL pwint_sl(stars,atoms,sym,zsl1,zsl2,&
                        volsli,volintsli,cell,slab%nmtsl(1,i),stars%kv3(1,j),x)
          stfunint(j,i) =  x*stars%nstr(j)
       ENDDO  ! over 3D stars
    ENDDO     ! over vacuum%layers
    !
    ! Here, I reordered the stuff to save memory
    !
    DO  n = 1,ne
       z_z(:) = CMPLX(0.0,0.0)
       q1 = 0.0
       IF (zmat%l_real) THEN
          DO  i = 1,lapw%nv(isp)
             q1 = q1 + zMat%z_r(i,n)*zMat%z_r(i,n)
          ENDDO
       ELSE
          DO  i = 1,lapw%nv(isp)
             q1 = q1 + REAL(zMat%z_c(i,n)*CONJG(zMat%z_c(i,n)))
          ENDDO
       ENDIF
       z_z(1) = q1/cell%omtil
       !
       !     ----> g.ne.0 stars
       !
       DO  i = 1,lapw%nv(isp)
          DO  j = 1,i-1
             ix1 = lapw%k1(j,isp) - lapw%k1(i,isp)
             iy1 = lapw%k2(j,isp) - lapw%k2(i,isp)
             iz1 = lapw%k3(j,isp) - lapw%k3(i,isp)
             IF (iabs(ix1).GT.stars%mx1) CYCLE
             IF (iabs(iy1).GT.stars%mx2) CYCLE
             IF (iabs(iz1).GT.stars%mx3) CYCLE
             ind = stars%ig(ix1,iy1,iz1)
             indp = stars%ig(-ix1,-iy1,-iz1)
             IF (ind.EQ.0 .OR. indp.EQ.0) CYCLE
             phase = stars%rgphs(ix1,iy1,iz1)/ (stars%nstr(ind)*cell%omtil)
             phasep = stars%rgphs(-ix1,-iy1,-iz1)/ (stars%nstr(indp)*cell%omtil)
             IF (zmat%l_real) THEN
                z_z(ind)  = z_z(ind)  + zMat%z_r(j,n)*zMat%z_r(i,n)*REAL(phase)
                z_z(indp) = z_z(indp) + zMat%z_r(i,n)*zMat%z_r(j,n)*REAL(phasep)
             ELSE
                z_z(ind) = z_z(ind) +zMat%z_c(j,n)*CONJG(zMat%z_c(i,n))*phase
                z_z(indp)= z_z(indp)+zMat%z_c(i,n)*CONJG(zMat%z_c(j,n))*phasep
             ENDIF
          ENDDO
       ENDDO
       ! ----> calculate a charge in the layer interstitial region of the film
       !
       DO i = 1,slab%nsl
          qi = 0.0
          DO j = 1,stars%ng3
             qi = qi + z_z(j)*stfunint(j,i)
          ENDDO
          slab%qintsl(i,n,ikpt,isp) = qi 
       ENDDO    ! over vacuum%layers         

    ENDDO ! over states

    DEALLOCATE ( stfunint, z_z ) 

  END SUBROUTINE q_int_sl
END MODULE m_qintsl
