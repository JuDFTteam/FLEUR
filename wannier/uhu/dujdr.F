      MODULE m_dujdr
      CONTAINS

      SUBROUTINE dujdr(jmtd,jri,rmsh,dx,ub,j2,b,l,lmaxd,dub)

      USE m_constants
      USE m_difcub

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: jri,jmtd,lmaxd
      REAL, INTENT(IN) :: rmsh(jmtd),dx
      REAL, INTENT(IN) :: ub(jmtd,2)       ! u(b2)
      REAL, INTENT(IN) :: j2(0:lmaxd,jmtd) ! j_l(b2*r)
      REAL, INTENT(IN) :: b                ! b2
      INTEGER, INTENT(IN) :: l           ! l of sph. Bessel j2

      REAL, INTENT(OUT) :: dub(jmtd,2)
      REAL :: xi,t(jri,2)
      INTEGER :: i,j

      ! derivatives d/dr for large and small component of q
      DO i=1,jri
       t(i,:) = ub(i,:) / rmsh(i) * j2(l,i)
      ENDDO

      DO j = 1, 2
         ! derivative at 1st point
         dub(1,j) = difcub( rmsh(1),t(1,j),rmsh(1) )

         ! derivative at 2nd...(jri-2)th point
         DO i = 2, jri-2
            dub(i,j) = difcub( rmsh(i-1),t(i-1,j),rmsh(i) )
         ENDDO

         ! derivative at last two points
         dub(jri-1,j) = difcub( rmsh(jri-3),t(jri-3,j),rmsh(jri-1) )
         dub(jri,j) = difcub( rmsh(jri-3),t(jri-3,j),rmsh(jri) )
      ENDDO

      DO i=1,jri
       dub(i,:) = dub(i,:)*rmsh(i)
      ENDDO

      ! complete d/dr (ub*j2) = ub'j2 + ub j2' with sph. Bessel func. j
      ! rule:  j'_{l}(ar) = a*j_{l-1}(ar) - (l+1)/r*j_{l}(ar)
c      IF(l.ne.0) THEN
c       DO i=1,jri
c        xi = rmsh(i)
c        dub(i,:) = dub(i,:) * j2(l,i) * xi
c     >           + ub(i,:) *( j2(l-1,i)*b - (l+1)/xi*j2(l,i) ) 
c       ENDDO
c      ELSE
c       DO i=1,jri
c        xi = rmsh(i)
c        dub(i,:) = dub(i,:) * j2(l,i) * xi
c     >           - ub(i,:) * j2(1,i) * b
c       ENDDO
c      ENDIF
      
      END SUBROUTINE dujdr
      END MODULE m_dujdr
