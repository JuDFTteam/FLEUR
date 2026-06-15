      MODULE m_clebsch
      CONTAINS
      REAL FUNCTION clebsch(aj,bj,am,bm,cj,cm)
******************************************************************
*     Program calculates Clebsch-Gordan coefficients                 *
*     See: Landau and Lifshitz, Vol.3                                *
*     cj,cm                                                     *
*     C                                                          *
*     aj,am,bj,bm                                               *
*     Written by A.Soldatov (IAE)                                    *
******************************************************************
      IMPLICIT NONE

      REAL, INTENT (IN) :: aj,bj,am,bm,cj,cm

      INTEGER n,k,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
      INTEGER iaj2,ibj2,icj2,im2,in2,ip2
      REAL    x,s,c,e,tol
      REAL    f(100)
      INTRINSIC sqrt,exp,min0,max0,nint,abs

      clebsch = 0.e0

      tol = 1.e-5
      iaj2 = NINT(2.e0*aj)
      ibj2 = NINT(2.e0*bj)
      icj2 = NINT(2.e0*cj)
      im2  = NINT(2.e0*am)
      in2  = NINT(2.e0*bm)
      ip2  = NINT(2.e0*cm)

      IF (ABS(2.e0*aj - REAL(iaj2)).GT.tol) RETURN
      IF (ABS(2.e0*bj - REAL(ibj2)).GT.tol) RETURN
      IF (ABS(2.e0*cj - REAL(icj2)).GT.tol) RETURN
      IF (ABS(2.e0*am - REAL(im2 )).GT.tol) RETURN
      IF (ABS(2.e0*bm - REAL(in2 )).GT.tol) RETURN
      IF (ABS(2.e0*cm - REAL(ip2 )).GT.tol) RETURN

      IF (im2 + in2 .NE. ip2) RETURN

      n = 100
      k =   0
      x = 2.0
      f(1)=0.e0 ; f(2)=0.e0

      IF ( k <= 0 ) THEN
         k=1
         DO i = 3,n
            f(i)=f(i-1)+log(x)
            x=x+1.e0
         ENDDO
      ENDIF

      i = am + bm - cm + .1e0
      IF (( i < 0 ).OR.( i > 0 )) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i1 = aj + bj - cj + 1.1e0
      IF ( i1 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i2 = aj - bj + cj + 1.1e0
      IF ( i2 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i3 = bj + cj - aj + 1.1e0
      IF ( i3 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      x  = aj + bj + cj + 2.1e0
      i4 = x
      i  = x + .6e0
      i  = i4 - i
      IF (( i < 0 ).OR.( i > 0 )) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      x = aj + am + 1.1e0
      i5 = x
      IF (i5 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i = x + .6e0
      i = i - I5
      IF (( i < 0 ).OR.( i > 0 )) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i6 = aj - am + 1.1e0
      IF (i6 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      x = bj + bm + 1.1e0
      i7=X
      IF (i7 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i = x + .6e0
      i = i - i7
      IF (( i < 0 ).OR.( i > 0 )) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i8 = bj - bm + 1.1e0
      IF (i8 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      x = cj + cm + 1.1e0
      i9 = x
      IF (i9 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i = x + .6e0
      i = i - i9
      IF (( i < 0 ).OR.( i > 0 )) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      i10 = cj - cm + 1.1e0
      IF (i10 <= 0 ) THEN
         clebsch = 0.e0
         RETURN
      ENDIF

      x = f(i1) + f(i2) + f(i3) - f(i4)
      i = i5 - i6
      IF ( i == 0 ) THEN
         i = i7 - i8
         IF ( i == 0 ) THEN
            i  = i4/2
            i5 = i4*0.5e0 + 0.6e0
            i = i - i5
            IF (( i < 0 ).OR.( i > 0 )) THEN
               clebsch = 0.e0
               RETURN
            ENDIF

            i6 = i5 - i6 + 1
            i7 = i5 - i8 + 1
            i8 = i5 - i10 + 1
            s = x*0.5e0 + f(i5) - f(i6) - f(i7) - f(i8)
            s = exp(s)
            i5 = i8/2
            i6 = i8*0.5e0 + 0.6e0
            i5 = i5 - I6
            IF ( i5 == 0 ) THEN
               s = 1.e0 - s - 1.e0
            ENDIF
            
            clebsch = s*SQRT( cj + cj + 1.e0 )
            return
         ENDIF
      ENDIF
      
      x = x + f(i5) + f(i6) + f(i7) + f(i8) + f(i9) + f(i10)
      x = x*0.5e0
      i10 = MIN0(i1,i6,i7)
      i2 = i1 - i5
      i3 = i1 - i8
      i9 = MAX0(0,i2,i3) + 1
      i1 = i1 + 1
      i6 = i6 + 1
      i7 = i7 + 1
      i8 = i9/2
      e  = 1.e0
      i5 = i9*0.5e0 + 0.6e0
      i8 = i8 - i5

      IF ( i8 == 0 ) THEN
         e  = -1.e0
      ENDIF
      S = 0.e0
      DO i = i9, i10
         c = x-f(i)-f(i1-i)-f(i6-i)-f(i7-i)-f(i-i2)-f(i-i3)
         s = s + e*exp(c)
         e = 1.e0 - e - 1.e0
      ENDDO
      clebsch = SQRT( cj + cj + 1.e0 )*s
      RETURN

      

      END FUNCTION clebsch

      SUBROUTINE clebsch_selftest(unit)
      IMPLICIT NONE
      INTEGER, INTENT(IN), OPTIONAL :: unit
      INTEGER out
      REAL c1,c2,c3,c4

      out = 6
      IF (PRESENT(unit)) out = unit

      c1 = clebsch(1.0,0.5,-1.0, 0.5,1.5,-0.5)
      c2 = clebsch(1.0,0.5, 0.0, 0.5,1.5, 0.5)
      c3 = clebsch(1.0,0.5, 1.0, 0.5,1.5, 1.5)
      c4 = clebsch(1.0,0.5, 0.0, 0.5,1.5, 1.5)

      WRITE(out,'(a)') 'Clebsch-Gordan self-test table:'
      WRITE(out,'(a,3f10.5,a,es18.10,a,es18.10)')
     + '<1,-1;1/2,+1/2|3/2,-1/2> ',1.0,-1.0,-0.5,
     + ' value=',c1,' expected=',1.0/SQRT(3.0)
      WRITE(out,'(a,3f10.5,a,es18.10,a,es18.10)')
     + '<1, 0;1/2,+1/2|3/2,+1/2> ',1.0,0.0,0.5,
     + ' value=',c2,' expected=',SQRT(2.0/3.0)
      WRITE(out,'(a,3f10.5,a,es18.10,a,es18.10)')
     + '<1,+1;1/2,+1/2|3/2,+3/2> ',1.0,1.0,1.5,
     + ' value=',c3,' expected=',1.0
      WRITE(out,'(a,3f10.5,a,es18.10,a,es18.10)')
     + '<1, 0;1/2,+1/2|3/2,+3/2> ',1.0,0.0,1.5,
     + ' value=',c4,' expected=',0.0
      END SUBROUTINE clebsch_selftest
      END MODULE m_clebsch
