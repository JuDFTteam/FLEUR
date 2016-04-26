      MODULE m_potl0
c     ******************************************************************
c     evaluate the xc-potential vxc for charge density and its
c     gradients,dens,... only for nonmagnetic.
c     ******************************************************************
      CONTAINS
      SUBROUTINE potl0(
     >                 mshd,jspd,jspins,icorr,msh,dx,rad,dens,
     <                 vxc)

      USE m_grdchlh
      USE m_mkgl0
      USE m_xcallg , ONLY : vxcallg

      IMPLICIT NONE
c     ..
      INTEGER, INTENT (IN) :: jspins,jspd,mshd,msh
      INTEGER, INTENT (IN) :: icorr
      REAL,    INTENT (IN) :: dx
      REAL,    INTENT (IN) :: rad(msh),dens(mshd,jspd)
      REAL,    INTENT (OUT):: vxc(mshd,jspd)

c     .. previously untyped names ..
      INTEGER,PARAMETER :: ndvgrd=6

      INTEGER i,ispin
      REAL, ALLOCATABLE :: drr(:,:),ddrr(:,:),agrt(:),agrd(:),agru(:)
      REAL, ALLOCATABLE :: g2rt(:),g2rd(:),g2ru(:),gggrt(:),gggrd(:)
      REAL, ALLOCATABLE :: gggru(:),grgrd(:),grgru(:),gzgr(:)
      
      REAL              :: vx(mshd,jspd)

      ALLOCATE ( drr(mshd,jspd),ddrr(mshd,jspd),grgru(mshd),gzgr(mshd) )
      ALLOCATE ( agrt(mshd),agrd(mshd),agru(mshd),g2rt(mshd),g2rd(mshd),
     +      g2ru(mshd),gggrt(mshd),gggrd(mshd),gggru(mshd),grgrd(mshd) )

      agrt(:) = 0.0 ; agru(:) = 0.0 ; agrd(:) = 0.0 ; grgrd(:) = 0.0
      g2rt(:) = 0.0 ; g2ru(:) = 0.0 ; g2rd(:) = 0.0 ; gzgr(:) = 0.0
      gggrt(:) = 0.0 ; gggru(:) = 0.0 ; gggrd(:) = 0.0 ; grgru(:) = 0.0 
!
!-->  evaluate gradients of dens.
!
      DO ispin = 1, jspins
        CALL grdchlh(
     >               1,1,msh,dx,rad,dens(1,ispin),ndvgrd,
     <               drr(1,ispin),ddrr(1,ispin))
      ENDDO

      CALL mkgl0(
     >           mshd,msh,jspd,jspins,rad,dens,drr,ddrr,
     <           agrt,agru,agrd,g2rt,g2ru,g2rd,
     <           gggrt,gggru,gggrd,grgru,grgrd,gzgr)
!
! --> calculate the potential.
!
      CALL vxcallg(
     >             icorr,.false.,jspins,mshd,msh,dens,
     +             agrt,agru,agrd,g2rt,g2ru,g2rd,gggrt,gggru,gggrd,
     +             gzgr,vx,vxc)

      DEALLOCATE ( drr,ddrr,grgru,gzgr,agrt,agrd,agru,g2rt,g2rd,
     +             g2ru,gggrt,gggrd,gggru,grgrd )
      END SUBROUTINE potl0
      END MODULE m_potl0
