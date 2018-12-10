c*************************************************c
c   compute radial integrals of sph. Hamiltonian  c
c   in scalar relativistic approximation between  c
c   products of radial function u and spherical   c
c   Bessel functions                              c
c                                                 c
c   integral = < uj1 | H_sra(lzb) | uj2 >         c
c*************************************************c
c                       J.-P. Hanke, Dec. 2015    c
c*************************************************c
      MODULE m_wann_uHu_radintsra5
      CONTAINS

      SUBROUTINE wann_uHu_radintsra5(jmtd,jri,rmsh,dx,e,vr,
     >                           uj1,uj2,duj1,duj2,lmaxd,lzb,
     >                           integral,irank)

      USE m_intgr, ONLY : intgr3
      USE m_constants
      IMPLICIT NONE

      REAL difcub
      EXTERNAL difcub
      
      INTEGER, INTENT(IN) :: jri,jmtd,lmaxd,irank
      REAL, INTENT(OUT) :: integral
      REAL, INTENT(IN) :: vr(jmtd),rmsh(jmtd),dx,e
      REAL, INTENT(IN) :: uj1(jmtd,2)      ! u(b1)*j
      REAL, INTENT(IN) :: uj2(jmtd,2)      ! u(b2)*j
      REAL, INTENT(IN) :: duj1(jmtd,2)     ! d/dr (u(b1)j)
      REAL, INTENT(IN) :: duj2(jmtd,2)     ! d/dr (u(b2)j)
      INTEGER, INTENT(IN) :: lzb         ! l of zentrifugal barrier

      REAL, ALLOCATABLE :: x(:)
      REAL :: c,c2,cin2,cin
      REAL :: ll,xi,vv,mm
      REAL :: sfp,facr,facl
      REAL, ALLOCATABLE :: xx(:,:),intt(:)
      INTEGER :: i,j,symopt

      c = c_light(1.)
      c2 = c*c
      cin = 1./c
      cin2 = cin*cin
      ll = lzb*(lzb+1)
      sfp = sqrt(4.0*pimach())
      symopt=0
      
      if(symopt.eq.0) then ! symmetrize
       facr = 0.5
       facl = 0.5
      elseif(symopt.lt.0)then ! apply H to left
       facr = 0.0
       facl = 1.0
      else ! apply H to right
       facr = 1.0
       facl = 0.0
      endif

      allocate( x(jri) )
      allocate( xx(jri,8),intt(8) )

      ! compute matrix elements of semi-relativistic
      ! Hamiltonian [Eq.(3.54) in PhD thesis of P.Kurz]
      DO i = 1, jri
         xi = rmsh(i)
         vv = vr(i) / xi !* sfp   ! no need to correct for Y_00
         mm = 1. + 0.5 * cin2 * ( e - vv )
         x(i) = 
     >        ! large-H-large
     >          uj1(i,1) * uj2(i,1) * ( 0.5 / mm * ll / xi / xi + vv )
     >        ! small-H-small
     >        + uj1(i,2) * uj2(i,2) * ( -2. * c2 + vv )
c      ! NEW VERSION BELOW
c     >        - c / xi * ( uj1(i,1)*uj2(i,2) + uj1(i,2)*uj2(i,1) )
c     >        + c / 2. * ( uj1(i,2)*duj2(i,1) + duj1(i,2)*uj2(i,1) )
c     >        - c / 2. * ( uj1(i,1)*duj2(i,2) + duj1(i,1)*uj2(i,2) )
      ! OLD VERSION BELOW
     >        ! large-H-small (symmetrized)
     >        - facr * c * uj1(i,1) * ( 2. * uj2(i,2) / xi + duj2(i,2) )
     >        - facl * c * uj2(i,1) * ( 2. * uj1(i,2) / xi + duj1(i,2) )
     >        ! small-H-large (symmetrized)
     >        + facr * c * uj1(i,2) * duj2(i,1)
     >        + facl * c * uj2(i,2) * duj1(i,1)

c         xx(i,:) = 0.0
c
c         xx(i,1)=uj1(i,1) * uj2(i,1) * ( 0.5 / mm * ll / xi / xi )
c         xx(i,2)=uj1(i,1) * uj2(i,1) * ( vv )
c
c         xx(i,3)=uj1(i,2) * uj2(i,2) * ( -2. * c2 )
c         xx(i,4)=uj1(i,2) * uj2(i,2) * ( vv )
c
c         xx(i,5)=-uj1(i,1) * c * ( 2.*uj2(i,2)/xi )
c         xx(i,6)=-uj1(i,1) * c * ( duj2(i,2) )
c
c         xx(i,7)= c * uj1(i,2) * duj2(i,1)
c
c         x(i) = xx(i,1)+xx(i,2)+xx(i,3)+xx(i,4)
c     >         +xx(i,5)+xx(i,6)+xx(i,7)+xx(i,8)

      ENDDO
      call intgr3(x,rmsh,dx,jri,integral)
c      do i=1,8
c       call intgr3(xx(:,i),rmsh,dx,jri,intt(i))
c      enddo
      !integral = integral / sfp

c      if(irank.eq.0 .and. abs(integral).gt.1e-10) then
c       do i=1,8
c        write(*,'(a,i1,a,f8.2)')'int',i,':',intt(i)/integral*100.
c       enddo
c       write(*,*)'integral:',integral
c      endif

      deallocate( x )
      deallocate( xx,intt )

      END SUBROUTINE wann_uHu_radintsra5
      END MODULE m_wann_uHu_radintsra5
