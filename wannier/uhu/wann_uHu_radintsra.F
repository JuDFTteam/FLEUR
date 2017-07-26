      MODULE m_wann_uHu_radintsra
      CONTAINS
      SUBROUTINE wann_uHu_radintsra(jmtd,jri,rmsh,dx,
     >                          epar,vr,f,g,l,expect)

      USE m_intgr, ONLY : intgr3
      USE m_constants
      IMPLICIT NONE

      REAL difcub
      EXTERNAL difcub
      
      REAL, INTENT(IN) :: f(jmtd,2),g(jmtd,2)
      REAL, INTENT(IN) :: vr(jmtd)
      REAL, INTENT(IN) :: epar
      REAL, INTENT(IN) :: rmsh(jmtd)
      REAL, INTENT(IN) :: dx
      REAL, INTENT(IN) :: expect
      INTEGER, INTENT(IN) :: jri
      INTEGER, INTENT(IN) :: jmtd
      INTEGER, INTENT(IN) :: l

      REAL, ALLOCATABLE :: x(:),dg(:,:),t(:,:),vv(:)
      REAL :: t11,t22,t12,t21,total,norm
      REAL :: mm,c,c2,cin2,cin
      REAL :: ll,xi,sfp
      INTEGER :: i,j

      c = c_light(1.)
      c2 = c*c
      cin = 1./c
      cin2 = cin*cin
      ll = l*(l+1)
      sfp = sqrt(4.0*pimach())

      allocate( x(jri), dg(jri,2), t(jri,2), vv(jri) )

      ! derivatives d/dr g for large and small component
      DO i=1,jri
       t(i,:) = g(i,:)/rmsh(i)
      ENDDO

      DO j = 1, 2
         ! derivative at 1st point
         dg(1,j) = difcub( rmsh(1),t(1,j),rmsh(1) )

         ! derivative at 2nd...(jri-2)th point
         DO i = 2, jri-2
            dg(i,j) = difcub( rmsh(i-1),t(i-1,j),rmsh(i) )
         ENDDO

         ! derivative at last two points
         dg(jri-1,j) = difcub( rmsh(jri-3),t(jri-3,j),rmsh(jri-1) )
         dg(jri,j) = difcub( rmsh(jri-3),t(jri-3,j),rmsh(jri) )
      ENDDO

      DO i=1,jri
       xi = rmsh(i)
       dg(i,:) = dg(i,:) * xi
       vv(i) = vr(i) / xi !* sfp
      ENDDO

      ! check normalization
      DO i = 1, jri
         x(i) = f(i,1)*g(i,1)+f(i,2)*g(i,2)
      ENDDO
      call intgr3(x,rmsh,dx,jri,norm)
      write(*,*)'norm:',norm


      ! compute matrix elements of semi-relativistic
      ! Hamiltonian [Eq.(3.54) in PhD thesis of P.Kurz]
      DO i = 1, jri
         mm = 1. + 0.5 * cin2 * ( epar - vv(i) )
         x(i) = f(i,1) * g(i,1)
     >        * ( 0.5 / mm * ll / rmsh(i) / rmsh(i) + vv(i) )
      ENDDO
      call intgr3(x,rmsh,dx,jri,t11) ! large-H-large

      DO i = 1, jri
         x(i) = f(i,2) * g(i,2) * ( -2. * c2 + vv(i) )
      ENDDO
      call intgr3(x,rmsh,dx,jri,t22) ! small-H-small

      DO i = 1, jri
         x(i) = f(i,1) * ( 2. * g(i,2) / rmsh(i) + dg(i,2) )
      ENDDO
      x = -c * x
      call intgr3(x,rmsh,dx,jri,t12) ! large-H-small

      DO i = 1, jri
         x(i) = f(i,2) * dg(i,1)
      ENDDO
      x = c * x
      call intgr3(x,rmsh,dx,jri,t21) ! small-H-large

      total = t11 + t22 + t12 + t21
      write(*,'(a,f16.12)')'expect:',expect
      write(*,'(a,f16.12)')'integr:',total
      if(abs(expect).gt.1e-12) then
       write(*,'(a,f7.3,a)')'differ:',abs(expect-total)/expect*100.,' %'
      endif

      deallocate( x, dg, t, vv )

      END SUBROUTINE wann_uHu_radintsra
      END MODULE m_wann_uHu_radintsra
