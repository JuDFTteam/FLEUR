      MODULE m_fitchk
      CONTAINS
      SUBROUTINE fitchk(f1,f2,av,rms,dmx)
!     ************************************************
!     compare functions f1 and f2
!     ************************************************
      IMPLICIT NONE
      REAL,INTENT(OUT):: av,dmx,rms
      REAL,INTENT(IN):: f1(:),f2(:)
!     .. Local Scalars ..
      REAL d
      INTEGER i

      av = 0.
      rms = 0.
      dmx = 0.
      DO  i = 1,SIZE(f1)
         av = av + f1(i)
         d = (f1(i)-f2(i))**2
         dmx = MAX(d,dmx)
         rms = rms + d
      ENDDO
      av = av/size(f1)
      IF (abs(av).LT.1.e-30) THEN
         rms = 0.
         dmx = 0.
         RETURN
      END IF
      rms = sqrt(rms/size(f1))/av*100.
      dmx = sqrt(dmx)/av*100.
      END SUBROUTINE 
      END
