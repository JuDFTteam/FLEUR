MODULE m_potl0
! ******************************************************************
! evaluate the xc-potential vxc for charge density and its
! gradients,dens,... only for nonmagnetic.
! ******************************************************************
CONTAINS
   SUBROUTINE potl0(xcpot,mshd,jspd,jspins,msh,dx,rad,dens, &
                    vxc)

      USE m_grdchlh
      USE m_mkgl0
      USE m_types
      IMPLICIT NONE

      CLASS(t_xcpot),intent(in)::xcpot
      INTEGER, INTENT (IN) :: jspins,jspd,mshd,msh
      REAL,    INTENT (IN) :: dx
      REAL,    INTENT (IN) :: rad(msh),dens(mshd,jspd)
      REAL,    INTENT (OUT):: vxc(mshd,jspd)

!     .. previously untyped names ..
      INTEGER,PARAMETER :: ndvgrd=6

      TYPE(t_gradients)::grad

      INTEGER i,ispin
      REAL, ALLOCATABLE :: drr(:,:),ddrr(:,:)

      REAL              :: vx(mshd,jspd)

      ALLOCATE ( drr(mshd,jspd),ddrr(mshd,jspd))
!
!-->  evaluate gradients of dens.
!
      CALL xcpot%alloc_gradients(msh,jspins,grad)
      DO ispin = 1, jspins
         CALL grdchlh(1,1,msh,dx,rad,dens(1,ispin),ndvgrd,&
                     drr(1,ispin),ddrr(1,ispin))
      ENDDO

      CALL mkgl0(mshd,msh,jspd,jspins,rad,dens,drr,ddrr,&
                 grad)
!
! --> calculate the potential.
!
      CALL xcpot%get_vxc(jspins, dens(:msh,:), vxc, vx, grad)

      DEALLOCATE ( drr,ddrr )
   END SUBROUTINE potl0
END MODULE m_potl0
