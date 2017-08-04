      MODULE m_d2fdz2cmplx
      CONTAINS

      SUBROUTINE d2fdz2cmplx(jmtd,jri,rmsh,dx,f,fac,d2f)

      IMPLICIT NONE

      REAL difcub
      EXTERNAL difcub
      
      COMPLEX, INTENT(INOUT) :: d2f(jmtd)
      COMPLEX, INTENT(IN)    :: fac(jmtd)

      REAL, INTENT(IN)    :: f(jmtd)
      REAL, INTENT(IN)    :: rmsh(jmtd)
      REAL, INTENT(IN)    :: dx

      INTEGER, INTENT(IN) :: jri
      INTEGER, INTENT(IN) :: jmtd

      REAL, ALLOCATABLE :: fr(:),fi(:),dfr(:),dfi(:)
      INTEGER :: i

      allocate( dfr(jri),dfi(jri),fr(jri),fi(jri) )
      DO i=1,jri
       fr(i) = f(i)*real(fac(i))
       fi(i) = f(i)*aimag(fac(i))
      ENDDO

      dfr(1) = difcub( rmsh(1),fr(1),rmsh(1) )
      dfi(1) = difcub( rmsh(1),fi(1),rmsh(1) )
      DO i = 2, jri-2
         dfr(i) = difcub( rmsh(i-1),fr(i-1),rmsh(i) )
         dfi(i) = difcub( rmsh(i-1),fi(i-1),rmsh(i) )
      ENDDO
      dfr(jri-1) = difcub( rmsh(jri-3),fr(jri-3),rmsh(jri-1) )
      dfi(jri-1) = difcub( rmsh(jri-3),fi(jri-3),rmsh(jri-1) )
      dfr(jri) = difcub( rmsh(jri-3),fr(jri-3),rmsh(jri) )
      dfi(jri) = difcub( rmsh(jri-3),fi(jri-3),rmsh(jri) )


      d2f(1) = cmplx( difcub( rmsh(1),dfr(1),rmsh(1) ),
     >                difcub( rmsh(1),dfi(1),rmsh(1) ) )
      DO i = 2, jri-2
         d2f(i) = cmplx( difcub( rmsh(i-1),dfr(i-1),rmsh(i) ),
     >                   difcub( rmsh(i-1),dfi(i-1),rmsh(i) ) )
      ENDDO
      d2f(jri-1) = cmplx( difcub( rmsh(jri-3),dfr(jri-3),rmsh(jri-1) ),
     >                    difcub( rmsh(jri-3),dfi(jri-3),rmsh(jri-1) ) )
      d2f(jri) = cmplx( difcub( rmsh(jri-3),dfr(jri-3),rmsh(jri) ),
     >                  difcub( rmsh(jri-3),dfi(jri-3),rmsh(jri) ) )

      deallocate( dfr,dfi,fr,fi )

c      d2f = cmplx(0.,0.)
c      d2f(1) = (f(3)-2*f(2)+f(1))/dx/dx
c      do i=2,jri-1
c       d2f(i) = (f(i+1)-2*f(i)+f(i-1))/dx/dx
c      enddo
c      d2f(jri) = (f(jri-2)-2*f(jri-1)+f(jri))/dx/dx 

      END SUBROUTINE d2fdz2cmplx
      END MODULE m_d2fdz2cmplx
