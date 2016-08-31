      MODULE m_qsf
      CONTAINS
      SUBROUTINE qsf(h,y,z,ndim,isave)
c     ..................................................................
c     subroutine qsf
c     purpose
c        to integrate an equidistant table of function values.
c     description of parameters
c         h     - the increment of arguement values.
c         y     - the input vector of function values.  y must be
c                 dimensioned ndim.
c         z     - if isave = 0,z is the final value of the integration.
c                 if isave not 0,z is the resulting vector of integrals
c                 for each y value. in this case it must be dimensioned
c                 ndim.
c         ndim  - the number of points in the table.
c         isave - determines whether a single value or multiple results
c                 will be returned
c                  isave =  0  - only the final value of the integra-
c                                tion will be returned.
c                  isave ne 0  - the results of integration at each
c                                point in the table will be returned.
c     remarks
c        no action in case ndim less than 3.
c     subroutines and function subprograms required
c        none
c     method
c        integration is done by means of simpsons rule together with
c        newtons 3/8 rule or a combination of these two rules. trunca-
c        tion error is of order h**5 if ndim greater than 3. if ndim = 3
c        truncation error is of order h**4.
c     reference
c        (1) f.b.hildebrand, introduction to numerical analysis,
c            mcgraw-hill, new york/toronto/london, 1956, pp.71-76.
c        (2) r. zurmuehl, praktische mathematic fuer ingenieure und
c            physiker, springer, berlin/goettingen/heidelberg, 1963,
c            pp. 214-221.
c     ..................................................................
C     .. Scalar Arguments ..

      IMPLICIT NONE
      REAL h
      INTEGER isave,ndim
C     ..
C     .. Array Arguments ..
      REAL y(ndim),z(*)
C     ..
C     .. Local Scalars ..
      REAL aux,aux1,aux2,ht,sum1,sum2,val
      INTEGER i
      LOGICAL save
C     ..
      save = isave .NE. 0
      ht = h / 3.0
      IF ( ndim == 5 ) THEN
        GOTO 140
      ELSEIF ( ndim < 5 ) THEN
        GOTO 130
      ENDIF
c     ndim is greater than 5. preparations of integration loop
      sum1 = y(2) + y(2)
      sum1 = sum1 + sum1
      sum1 = ht* (y(1)+sum1+y(3))
      aux1 = y(4) + y(4)
      aux1 = aux1 + aux1
      aux1 = sum1 + ht* (y(3)+aux1+y(5))
      aux2 = ht* (y(1)+3.875e0* (y(2)+y(5))+2.625e0* (y(3)+y(4))+y(6))
      IF (.NOT.save) GO TO 30
   20 z(1) = 0.e0
      sum2 = y(5) + y(5)
      sum2 = sum2 + sum2
      sum2 = aux2 - ht* (y(4)+sum2+y(6))
      aux = y(3) + y(3)
      aux = aux + aux
      z(2) = sum2 - ht* (y(2)+aux+y(4))
      z(3) = sum1
      z(4) = sum2
   30 CONTINUE
      IF ( ndim <= 6 ) THEN
         GOTO 70
      ENDIF
c     integration loop
   40 DO 60 i = 7,ndim,2
         sum1 = aux1
         sum2 = aux2
         aux1 = y(i-1) + y(i-1)
         aux1 = aux1 + aux1
         aux1 = sum1 + ht* (y(i-2)+aux1+y(i))
         IF (save) z(i-2) = sum1
         IF ( i >= ndim ) THEN
           GOTO 100
         ENDIF
         aux2 = y(i) + y(i)
         aux2 = aux2 + aux2
         aux2 = sum2 + ht* (y(i-1)+aux2+y(i+1))
         IF (save) z(i-1) = sum2
   60 CONTINUE
   70 IF (.NOT.save) GO TO 90
   80 z(ndim-1) = aux1
      z(ndim) = aux2
      RETURN
   90 z(1) = aux2
      RETURN
  100 IF (.NOT.save) GO TO 120
  110 z(ndim-1) = sum2
      z(ndim) = aux1
      RETURN
  120 z(1) = aux1
      RETURN
c     end of integration loop
  130 CONTINUE
      IF ( ndim == 3 ) THEN
        GOTO 210
      ELSEIF ( ndim < 3 ) THEN
        GOTO 230
      ENDIF
c     ndim is equal to 4 or 5
  140 sum2 = 1.125e0*ht* (y(1)+y(2)+y(2)+y(2)+y(3)+y(3)+y(3)+y(4))
      sum1 = y(2) + y(2)
      sum1 = sum1 + sum1
      sum1 = ht* (y(1)+sum1+y(3))
      IF (.NOT.save) GO TO 160
  150 z(1) = 0.e0
      aux1 = y(3) + y(3)
      aux1 = aux1 + aux1
      z(2) = sum2 - ht* (y(2)+aux1+y(4))
  160 CONTINUE
      IF ( ndim < 5 ) THEN
        GOTO 190
      ENDIF
  170 aux1 = y(4) + y(4)
      aux1 = aux1 + aux1
      val = sum1 + ht* (y(3)+aux1+y(5))
      IF (save) GO TO 180
      z(1) = val
      RETURN
  180 z(5) = val
      GO TO 200
  190 IF (save) GO TO 200
      z(1) = sum2
      RETURN
  200 z(3) = sum1
      z(4) = sum2
      RETURN
c     ndim is equal to 3
  210 sum2 = y(2) + y(2)
      sum2 = sum2 + sum2
      val = ht* (y(1)+sum2+y(3))
      IF (save) GO TO 220
      z(1) = val
      RETURN
  220 z(1) = 0.e0
      z(2) = ht* (1.25e0*y(1)+y(2)+y(2)-.25e0*y(3))
      z(3) = val
  230 RETURN
      END SUBROUTINE qsf
      END MODULE m_qsf
