      MODULE m_d2fdz2
      CONTAINS

      SUBROUTINE d2fdz2(jmtd,jri,rmsh,dx,f,fac,d2f)

      USE m_difcub

      IMPLICIT NONE
      
      REAL, INTENT(INOUT) :: d2f(jmtd)
      REAL, INTENT(IN)    :: fac(jmtd)
      REAL, INTENT(IN)    :: f(jmtd)
      REAL, INTENT(IN)    :: rmsh(jmtd)
      REAL, INTENT(IN)    :: dx

      INTEGER, INTENT(IN) :: jri
      INTEGER, INTENT(IN) :: jmtd

      REAL, ALLOCATABLE :: fr(:),dfr(:)
      INTEGER :: i

      allocate( dfr(jri),fr(jri) )
      DO i=1,jri
       fr(i) = f(i)*fac(i)
      ENDDO

      dfr(1) = difcub( rmsh(1),fr(1),rmsh(1) )
      DO i = 2, jri-2
         dfr(i) = difcub( rmsh(i-1),fr(i-1),rmsh(i) )
      ENDDO
      dfr(jri-1) = difcub( rmsh(jri-3),fr(jri-3),rmsh(jri-1) )
      dfr(jri) = difcub( rmsh(jri-3),fr(jri-3),rmsh(jri) )


      d2f(1) = difcub( rmsh(1),dfr(1),rmsh(1) )
      DO i = 2, jri-2
         d2f(i) = difcub( rmsh(i-1),dfr(i-1),rmsh(i) )
      ENDDO
      d2f(jri-1) = difcub( rmsh(jri-3),dfr(jri-3),rmsh(jri-1) )
      d2f(jri) = difcub( rmsh(jri-3),dfr(jri-3),rmsh(jri) )


      deallocate( dfr,fr )

      END SUBROUTINE d2fdz2
      END MODULE m_d2fdz2
