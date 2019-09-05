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
      MODULE m_wann_uHu_radintsra4
      CONTAINS

      SUBROUTINE wann_uHu_radintsra4(jmtd,jri,rmsh,dx,e,vr,vm,
     >                           uj1,uj2,duj2,lmaxd,lzb,
     >                           integral)

      USE m_intgr, ONLY : intgr3
      USE m_constants
      IMPLICIT NONE

      REAL difcub
      EXTERNAL difcub
      
      INTEGER, INTENT(IN) :: jri,jmtd,lmaxd
      REAL, INTENT(OUT) :: integral
      REAL, INTENT(IN) :: vr(jmtd),vm(jmtd),rmsh(jmtd),dx,e
      REAL, INTENT(IN) :: uj1(jmtd,2)      ! u(b1)*j
      REAL, INTENT(IN) :: uj2(jmtd,2)      ! u(b2)*j
      REAL, INTENT(IN) :: duj2(jmtd,2)     ! d/dr (u(b2)j)
      INTEGER, INTENT(IN) :: lzb         ! l of zentrifugal barrier

      REAL, ALLOCATABLE :: x(:)
      REAL :: c,c2,cin2,cin
      REAL :: ll,xi,vv,mm,vv2
      REAL :: sfp
      INTEGER :: i,j

      c = c_light(1.)
      c2 = c*c
      cin = 1./c
      cin2 = cin*cin
      ll = lzb*(lzb+1)
      sfp = sqrt(4.0*pimach())

      allocate( x(jri) )

      ! compute matrix elements of semi-relativistic
      ! Hamiltonian [Eq.(3.54) in PhD thesis of P.Kurz]
      DO i = 1, jri
         xi = rmsh(i)
         vv = vr(i) / xi !* sfp
         vv2= vm(i) / xi !* sfp
         mm = 1. + 0.5 * cin2 * ( e - vv2 )
         x(i) = 
     >        ! large-H-large
     >          uj1(i,1) * uj2(i,1) * ( 0.5 / mm * ll / xi / xi + vv )
     >        ! small-H-small
     >        + uj1(i,2) * uj2(i,2) * ( -2. * c2 + vv )
     >        ! large-H-small (not symmetrized)
     >        - c * uj1(i,1) * (2. * uj2(i,2) / xi + duj2(i,2) )
     >        ! small-H-large (not symmetrized)
     >        + c * uj1(i,2) * duj2(i,1)
      ENDDO
      call intgr3(x,rmsh,dx,jri,integral)

      deallocate( x )

      END SUBROUTINE wann_uHu_radintsra4
      END MODULE m_wann_uHu_radintsra4
