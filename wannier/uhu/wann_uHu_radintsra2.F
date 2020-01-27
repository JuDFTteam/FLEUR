      MODULE m_wann_uHu_radintsra2
      CONTAINS

      ! < p | H(l) | q >
      SUBROUTINE wann_uHu_radintsra2(jmtd,jri,rmsh,dx,
     >                           e,vr,
     >                           p,q,l,
     >                           integral)

      USE m_intgr, ONLY : intgr3
      USE m_constants
      USE m_difcub

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: jmtd
      REAL, INTENT(IN) :: p(jmtd,2)
      REAL, INTENT(IN) :: q(jmtd,2)
      REAL, INTENT(IN) :: vr(jmtd)
      REAL, INTENT(IN) :: e
      REAL, INTENT(IN) :: rmsh(jmtd)
      REAL, INTENT(IN) :: dx
      INTEGER, INTENT(IN) :: jri
      INTEGER, INTENT(IN) :: l
      REAL, INTENT(OUT) :: integral

      REAL, ALLOCATABLE :: x(:),dq(:,:),t(:,:)
      REAL :: c,c2,cin2,cin
      REAL :: ll,xi,vv,mm,sfp
      INTEGER :: i,j

      c = c_light(1.)
      c2 = c*c
      cin = 1./c
      cin2 = cin*cin
      ll = l*(l+1)
      sfp = sqrt(4.0*pimach())

      allocate( x(jri), dq(jri,2), t(jri,2) )

      DO i=1,jri
       t(i,:) = q(i,:) / rmsh(i)
      ENDDO

      ! derivatives d/dr for large and small component
      DO j = 1, 2
         ! derivative at 1st point
         dq(1,j) = difcub( rmsh(1),t(1,j),rmsh(1) )

         ! derivative at 2nd...(jri-2)th point
         DO i = 2, jri-2
            dq(i,j) = difcub( rmsh(i-1),t(i-1,j),rmsh(i) )
         ENDDO

         ! derivative at last two points
         dq(jri-1,j) = difcub( rmsh(jri-3),t(jri-3,j),rmsh(jri-1) )
         dq(jri,j) = difcub( rmsh(jri-3),t(jri-3,j),rmsh(jri) )
      ENDDO

      DO i=1,jri
       dq(i,:) = dq(i,:)*rmsh(i)
      ENDDO

      ! compute matrix elements of semi-relativistic
      ! Hamiltonian [Eq.(3.54) in PhD thesis of P.Kurz]
      DO i = 1, jri
         xi = rmsh(i)
         vv = vr(i) / xi !* sfp
         mm = 1. + 0.5 * cin2 * ( e - vv )
         x(i) = 
     >        ! large-H-large
     >          p(i,1) * q(i,1) * ( 0.5 / mm * ll / xi / xi + vv )
     >        ! small-H-small
     >        + p(i,2) * q(i,2) * ( -2. * c2 + vv )
     >        ! large-H-small (not symmetrized)
     >        - c * p(i,1) * (2. * q(i,2) / xi + dq(i,2) )
     >        ! small-H-large (not symmetrized)
     >        + c * p(i,2) * dq(i,1)
      ENDDO
      call intgr3(x,rmsh,dx,jri,integral)

      deallocate( x, dq, t )

      END SUBROUTINE wann_uHu_radintsra2
      END MODULE m_wann_uHu_radintsra2
