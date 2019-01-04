MODULE m_grdchlh
      use m_juDFT
!     -----------------------------------------------------------------
!     input: iexpm=1: exponential mesh. otherwise dx interval mesh.
!            ro: charge or quantity to be derivated.
!     evaluates d(ro)/dr,d{d(ro)/dr}/dr.
!     drr=d(ro)/dr, ddrr=d(drr)/dr.
!     coded by t.asada. 1996.
!     ------------------------------------------------------------------
      CONTAINS

      SUBROUTINE grdchlh(iexpm,ist,ied,dx,rv,ro,ndvgrd, drr,ddrr)

      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: iexpm,ist,ied,ndvgrd
      REAL,    INTENT (IN)  :: dx
      REAL,    INTENT (IN)  :: ro(ied),rv(ied)
      REAL,    INTENT (OUT) :: drr(ied),ddrr(ied)

      REAL drx,drx0,drx1,drx2,drx3,drxx,drxx0,drxx1,drxx2,drxx3
      INTEGER i,i1,i2,i3,i4,i5,i6,j,nred

      DO i = ist,ied
          drr(i) = 0.0
          ddrr(i) = 0.0
      ENDDO

      IF (ied-ist < 3) RETURN


      IF (ndvgrd < 3 .OR. ndvgrd.GT.6) THEN
          WRITE (16,fmt=126) ndvgrd
          CALL juDFT_error("ndvgrd<3 .or. ndvgrd>6",calledby="grdchlh")
      ENDIF
  126 FORMAT (/,' ndvgrd should be ge.4 .or. le.6. ndvgrd=',i3)

      i1 = ist
      i2 = ist + 1
      i3 = ist + 2
      i4 = ist + 3
      i5 = ist + 4
      i6 = ist + 5

      IF (ndvgrd==3) THEN

          drx1  = f131(ro(i1),ro(i2),ro(i3),dx)
          drxx1 = f231(ro(i1),ro(i2),ro(i3),dx)

      ELSEIF (ndvgrd==4) THEN

          drx1  = f141(ro(i1),ro(i2),ro(i3),ro(i4),dx)
          drxx1 = f241(ro(i1),ro(i2),ro(i3),ro(i4),dx)
          drx2  = f142(ro(i1),ro(i2),ro(i3),ro(i4),dx)
          drxx2 = f242(ro(i1),ro(i2),ro(i3),ro(i4),dx)

      ELSEIF (ndvgrd==5) THEN

          drx1  = f151(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
          drxx1 = f251(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
          drx2  = f152(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)
          drxx2 = f252(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),dx)

      ELSEIF (ndvgrd==6) THEN

          drx1  = f161(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drxx1 = f261(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drx2  = f162(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drxx2 = f262(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drx3  = f163(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)
          drxx3 = f263(ro(i1),ro(i2),ro(i3),ro(i4),ro(i5),ro(i6),dx)

      ENDIF

      IF (iexpm==1) THEN
          drr(i1)  = drx1/rv(i1)
          ddrr(i1) = (drxx1-drx1)/rv(i1)**2
      ELSE
          drr(i1)  = drx1
          ddrr(i1) = drxx1
      ENDIF

      IF (ndvgrd.GT.3) THEN

          IF (iexpm==1) THEN
              drr(i2) = drx2/rv(i2)
              ddrr(i2) = (drxx2-drx2)/rv(i2)**2
          ELSE
              drr(i2) = drx2
              ddrr(i2) = drxx2
          ENDIF

          IF (ndvgrd==6) THEN
              IF (iexpm==1) THEN
                  drr(i3)  = drx3/rv(i3)
                  ddrr(i3) = (drxx3-drx3)/rv(i3)**2
              ELSE
                  drr(i3)  = drx3
                  ddrr(i3) = drxx3
              ENDIF
          ENDIF

      ENDIF

      nred = REAL(ndvgrd)/2 + 0.1

      IF (ied-nred.LE.ist) THEN
          WRITE(16,fmt='(/'' ied-nred < ist. ied,nred,ist='',3i4)') ied,nred,ist
           CALL juDFT_error("ied-nred.le.ist",calledby="grdchlh")
      ENDIF

      DO j = nred + ist,ied - nred

          IF (ndvgrd==3) THEN

              drx  = f132(ro(j-1),ro(j),ro(j+1),dx)
              drxx = f232(ro(j-1),ro(j),ro(j+1),dx)

          ELSEIF (ndvgrd==4) THEN

              drx  = f142(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
              drxx = f242(ro(j-1),ro(j),ro(j+1),ro(j+2),dx)

          ELSEIF (ndvgrd==5) THEN

              drx  = f153(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)
              drxx = f253(ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2),dx)

          ELSEIF (ndvgrd==6) THEN

              drx  = f164(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2), dx)
              drxx = f264(ro(j-3),ro(j-2),ro(j-1),ro(j),ro(j+1),ro(j+2), dx)

          ENDIF

          IF (iexpm==1) THEN
              drr(j)  = drx/rv(j)
              ddrr(j) = (drxx-drx)/rv(j)**2
          ELSE
              drr(j)  = drx
              ddrr(j) = drxx
          ENDIF


      ENDDO ! j
      
      IF (ndvgrd==3) THEN

          drx0  = f133(ro(ied-2),ro(ied-1),ro(ied),dx)
          drxx0 = f233(ro(ied-2),ro(ied-1),ro(ied),dx)

      ELSEIF (ndvgrd==4) THEN

          drx1  = f143(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)
          drxx1 = f243(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)
          drx0  = f144(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)
          drxx0 = f244(ro(ied-3),ro(ied-2),ro(ied-1),ro(ied),dx)

      ELSEIF (ndvgrd==5) THEN

          drx1  = f154(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied), dx)
          drxx1 = f254(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied), dx)
          drx0  = f155(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied), dx)
          drxx0 = f255(ro(ied-4),ro(ied-3),ro(ied-2),ro(ied-1),ro(ied), dx)

      ELSEIF (ndvgrd==6) THEN

          drx2  = f164(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2), ro(ied-1),ro(ied),dx)
          drxx2 = f264(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2), ro(ied-1),ro(ied),dx)

          drx1  = f165(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2), ro(ied-1),ro(ied),dx)
          drxx1 = f265(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2), ro(ied-1),ro(ied),dx)

          drx0  = f166(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2), ro(ied-1),ro(ied),dx)
          drxx0 = f266(ro(ied-5),ro(ied-4),ro(ied-3),ro(ied-2), ro(ied-1),ro(ied),dx)


      ENDIF

      IF (ndvgrd.GT.3) THEN

          IF (ndvgrd==6) THEN
              IF (iexpm==1) THEN
                  drr(ied-2)  = drx2/rv(ied-2)
                  ddrr(ied-2) = (drxx2-drx2)/rv(ied-2)**2
              ELSE
                  drr(ied-2)  = drx2
                  ddrr(ied-2) = drxx2
              ENDIF
          ENDIF

          IF (iexpm==1) THEN
              drr(ied-1)  = drx1/rv(ied-1)
              ddrr(ied-1) = (drxx1-drx1)/rv(ied-1)**2
          ELSE
              drr(ied-1)  = drx1
              ddrr(ied-1) = drxx1
          ENDIF

      ENDIF

      IF (iexpm==1) THEN
          drr(ied)  = drx0/rv(ied)
          ddrr(ied) = (drxx0-drx0)/rv(ied)**2
      ELSE
          drr(ied)  = drx0
          ddrr(ied) = drxx0
      ENDIF

      END SUBROUTINE grdchlh
!--------------------------------------------------------------------
! Functions: 3 point formula :
!
      REAL FUNCTION f131(f0,f1,f2,d)
        REAL f0,f1,f2,d
        f131 = (-3*f0+4*f1-f2)/ (2*d)
      END FUNCTION f131
      REAL FUNCTION f132(g1,f0,f1,d)
        REAL g1,f0,f1,d
        f132 = (-1*g1-0*f0+f1)/ (2*d)
      END FUNCTION f132
      REAL FUNCTION f133(g2,g1,f0,d)
        REAL g2,g1,f0,d
        f133 = (g2-4*g1+3*f0)/ (2*d)
      END FUNCTION f133
!
!.....four point formula for the 1st deriv.
      REAL FUNCTION f141(f0,f1,f2,f3,d)
        REAL f0,f1,f2,f3,d
        f141 = (-11*f0+18*f1-9*f2+2*f3)/ (6*d)
      END FUNCTION f141
      REAL FUNCTION f142(g1,f0,f1,f2,d)
        REAL g1,f0,f1,f2,d
        f142 = (-2*g1-3*f0+6*f1-f2)/ (6*d)
      END FUNCTION f142
      REAL FUNCTION f143(g2,g1,f0,f1,d)
        REAL g2,g1,f0,f1,d
        f143 = (g2-6*g1+3*f0+2*f1)/ (6*d)
      END FUNCTION f143
      REAL FUNCTION f144(g3,g2,g1,f0,d)
        REAL g3,g2,g1,f0,d
        f144 = (-2*g3+9*g2-18*g1+11*f0)/ (6*d)
      END FUNCTION f144
!
!.....five point formula for the 1st deriv.
      REAL FUNCTION f151(f0,f1,f2,f3,f4,d)
        REAL f0,f1,f2,f3,f4,d 
        f151 = (-50*f0+96*f1-72*f2+32*f3-6*f4)/ (24*d)
      END FUNCTION f151
      REAL FUNCTION f152(g1,f0,f1,f2,f3,d)
        REAL g1,f0,f1,f2,f3,d
        f152 = (-6*g1-20*f0+36*f1-12*f2+2*f3)/ (24*d)
      END FUNCTION f152
      REAL FUNCTION f153(g2,g1,f0,f1,f2,d)
        REAL g2,g1,f0,f1,f2,d
        f153 = (2*g2-16*g1-0*f0+16*f1-2*f2)/ (24*d)
      END FUNCTION f153
      REAL FUNCTION f154(g3,g2,g1,f0,f1,d)
        REAL g3,g2,g1,f0,f1,d
        f154 = (-2*g3+12*g2-36*g1+20*f0+6*f1)/ (24*d)
      END FUNCTION f154
      REAL FUNCTION f155(g4,g3,g2,g1,f0,d)
        REAL g4,g3,g2,g1,f0,d
        f155 = (6*g4-32*g3+72*g2-96*g1+50*f0)/ (24*d)
      END FUNCTION f155
!
!.....six point formula for the 1st deriv.
      REAL FUNCTION f161(f0,f1,f2,f3,f4,f5,d)
        REAL f0,f1,f2,f3,f4,f5,d
        f161 = (-274*f0+600*f1-600*f2+400*f3-150*f4+24*f5)/ (120*d)
      END FUNCTION f161
      REAL FUNCTION f162(g1,f0,f1,f2,f3,f4,d)
        REAL g1,f0,f1,f2,f3,f4,d
        f162 = (-24*g1-130*f0+240*f1-120*f2+40*f3-6*f4)/ (120*d)
      END FUNCTION f162
      REAL FUNCTION f163(g2,g1,f0,f1,f2,f3,d)
        REAL g2,g1,f0,f1,f2,f3,d
        f163 = (6*g2-60*g1-40*f0+120*f1-30*f2+4*f3)/ (120*d)
      END FUNCTION f163
      REAL FUNCTION f164(g3,g2,g1,f0,f1,f2,d) 
        REAL g3,g2,g1,f0,f1,f2,d
        f164 = (-4*g3+30*g2-120*g1+40*f0+60*f1-6*f2)/ (120*d)
      END FUNCTION f164
      REAL FUNCTION f165(g4,g3,g2,g1,f0,f1,d)
        REAL g4,g3,g2,g1,f0,f1,d
        f165 = (6*g4-40*g3+120*g2-240*g1+130*f0+24*f1)/ (120*d)
      END FUNCTION f165
      REAL FUNCTION f166(g5,g4,g3,g2,g1,f0,d)
        REAL g5,g4,g3,g2,g1,f0,d
        f166 = (-24*g5+150*g4-400*g3+600*g2-600*g1+274*f0)/ (120*d)
      END FUNCTION f166
!
!.....three point formula for the 2nd deriv.
      REAL FUNCTION f231(f0,f1,f2,d)
        REAL f0,f1,f2,d
        f231 = (f0-2*f1+f2)/ (d*d)
      END FUNCTION f231
      REAL FUNCTION f232(g1,f0,f1,d)
        REAL g1,f0,f1,d
        f232 = (g1-2*f0+f1)/ (d*d)
      END FUNCTION f232
      REAL FUNCTION f233(g2,g1,f0,d)
        REAL g2,g1,f0,d
        f233 = (g2-2*g1+f0)/ (d*d)
      END FUNCTION f233
!
!.....four point formula for the 2nd deriv.
      REAL FUNCTION f241(f0,f1,f2,f3,d)
        REAL f0,f1,f2,f3,d
        f241 = (6*f0-15*f1+12*f2-3*f3)/ (3*d*d)
      END FUNCTION f241
      REAL FUNCTION f242(g1,f0,f1,f2,d)
        REAL g1,f0,f1,f2,d
        f242 = (3*g1-6*f0+3*f1+0*f2)/ (3*d*d)
      END FUNCTION f242
      REAL FUNCTION f243(g2,g1,f0,f1,d)
        REAL g2,g1,f0,f1,d
        f243 = (0*g2+3*g1-6*f0+3*f1)/ (3*d*d)
      END FUNCTION f243
      REAL FUNCTION f244(g3,g2,g1,f0,d)
        REAL g3,g2,g1,f0,d
        f244 = (-3*g3+2*g2+15*g1+6*f0)/ (3*d*d)
      END FUNCTION f244
!
!.....five point formula for the 2nd deriv.
      REAL FUNCTION f251(f0,f1,f2,f3,f4,d)
        REAL f0,f1,f2,f3,f4,d
        f251 = (35*f0-104*f1+114*f2-56*f3+11*f4)/ (12*d*d)
      END FUNCTION f251
      REAL FUNCTION f252(g1,f0,f1,f2,f3,d)
        REAL g1,f0,f1,f2,f3,d
        f252 = (11*g1-20*f0+6*f1+4*f2-f3)/ (12*d*d)
      END FUNCTION f252
      REAL FUNCTION f253(g2,g1,f0,f1,f2,d)
        REAL g2,g1,f0,f1,f2,d
        f253 = (-g2+16*g1-30*f0+16*f1-f2)/ (12*d*d)
      END FUNCTION f253
      REAL FUNCTION f254(g3,g2,g1,f0,f1,d)
        REAL g3,g2,g1,f0,f1,d
        f254 = (-g3+4*g2+6*g1-20*f0+11*f1)/ (12*d*d)
      END FUNCTION f254
      REAL FUNCTION f255(g4,g3,g2,g1,f0,d)
        REAL g4,g3,g2,g1,f0,d
        f255 = (11*g4-56*g3+114*g2-104*g1+35*f0)/ (12*d*d)
      END FUNCTION f255
!
!.....six point formula for the 2nd deriv.
      REAL FUNCTION f261(f0,f1,f2,f3,f4,f5,d)
        REAL f0,f1,f2,f3,f4,f5,d
        f261 = (225*f0-770*f1+1070*f2-780*f3+305*f4-50*f5)/ (60*d*d)
      END FUNCTION f261
      REAL FUNCTION f262(g1,f0,f1,f2,f3,f4,d)
        REAL g1,f0,f1,f2,f3,f4,d
        f262 = (50*g1-75*f0-20*f1+70*f2-30*f3+5*f4)/ (60*d*d)
      END FUNCTION f262
      REAL FUNCTION f263(g2,g1,f0,f1,f2,f3,d) 
        REAL g2,g1,f0,f1,f2,f3,d
        f263 = (-5*g2+80*g1-150*f0+80*f1-5*f2+0*f3)/ (60*d*d)
      END FUNCTION f263
      REAL FUNCTION f264(g3,g2,g1,f0,f1,f2,d)
        REAL g3,g2,g1,f0,f1,f2,d
        f264 = (0*g3-5*g2+80*g1-150*f0+80*f1-5*f2)/ (60*d*d)
      END FUNCTION f264
      REAL FUNCTION f265(g4,g3,g2,g1,f0,f1,d) 
        REAL g4,g3,g2,g1,f0,f1,d
        f265 = (5*g4-30*g3+70*g2-20*g1-75*f0+50*f1)/ (60*d*d)
      END FUNCTION f265
      REAL FUNCTION f266(g5,g4,g3,g2,g1,f0,d) 
        REAL g5,g4,g3,g2,g1,f0,d
        f266 = (-50*g5+305*g4-780*g3+1070*g2-770*g1+225*f0)/ (60*d*d)
      END FUNCTION f266
!--------------------------------------------------------------------
      END MODULE m_grdchlh
