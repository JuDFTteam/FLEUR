MODULE m_potl0
! ******************************************************************
! evaluate the xc-potential vxc for charge density and its
! gradients,dens,... only for nonmagnetic.
! ******************************************************************
CONTAINS
   SUBROUTINE potl0(xcpot,jspins,dx,rad,dens, &
                    vxc)

      USE m_grdchlh
      USE m_mkgl0
      USE m_types
      IMPLICIT NONE

      CLASS(t_xcpot),intent(in)::xcpot
      INTEGER, INTENT (IN) :: jspins
      REAL,    INTENT (IN) :: dx
      REAL,    INTENT (IN) :: rad(:),dens(:,:)
      REAL,    INTENT (OUT):: vxc(:,:)

!     .. previously untyped names ..

      TYPE(t_gradients)::grad

      INTEGER i,ispin,msh
      REAL, ALLOCATABLE :: drr(:,:),ddrr(:,:)

      REAL              :: vx(size(vxc,1),jspins)

      msh = size(rad)
      ALLOCATE ( drr(msh,jspins),ddrr(msh,jspins))
!
!-->  evaluate gradients of dens.
!
      CALL xcpot%alloc_gradients(msh,jspins,grad)
      DO ispin = 1, jspins
         CALL grdchlh(dx,dens(1:msh,ispin),&
                     drr(:,ispin),ddrr(:,ispin),rad)
      ENDDO

      CALL mkgl0(jspins,rad,dens,drr,ddrr,&
                 grad)
!
! --> calculate the potential.
!
      CALL xcpot%get_vxc(jspins, dens(:msh,:), vxc, vx, grad)

      DEALLOCATE ( drr,ddrr )
   END SUBROUTINE potl0
END MODULE m_potl0
