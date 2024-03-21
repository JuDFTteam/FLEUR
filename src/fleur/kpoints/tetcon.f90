MODULE m_tetcon
   use m_juDFT
   use m_constants

! This subroutine constructs the tetrahedra for the
! Brillouin zone integration

CONTAINS
   SUBROUTINE tetcon( &
      nkpt, ndiv3, &
      omega, kvc, nsym, &
      nt, voltet, ntetra)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ndiv3, nkpt, nsym
      REAL, INTENT(IN) :: omega
      REAL, INTENT(IN) :: kvc(3, nkpt)
      INTEGER, INTENT(OUT) :: nt
      INTEGER, INTENT(OUT) :: ntetra(4, ndiv3)
      REAL, INTENT(OUT) :: voltet(ndiv3)

      INTEGER :: i, it, j, ic, icom, ik, nkp, jj, nnt, nsid, nkq
      INTEGER :: i2, i3, i4, nk2, k, l, kk, is1, is2, is3, js1, js2, js3
      INTEGER :: nside(4), nside2(4), nkadd(nkpt)
      REAL :: dnorm, vav, vsq, sav, ssq, dm, dist, dl, vect
      REAL :: vol, vt, length, minlen, eps, eps1
      REAL :: vmin, vmax, smin, smax, volnew, sum, sss
      REAL :: ddist(nkpt), kcorn(3, 4), norm(3), d(3, 16), dnm(3), cn(3)

!----->  Intrinsic Functions

      INTRINSIC abs, max, min, sqrt

!*************************************************************************

!     All commons have been removed from the historical program and re-
!     placed by direct subroutine calls.  Those statements immediately
!     below preceded by stars your removed prior to my changes

!                                          June 9, 1992
!                                                   Fred Hutson

!*************************************************************************

      data eps/1e-9/
      data eps1/1e-5/

! CONSTRUCT THE FIRST TETRAHEDRON

      ntetra(1, 1) = 1

! FIND K-POINT NEAREST TO 1

      dm = 10.0**9
      i2 = 0
      do i = 2, nkpt
         dist = 0.0
         do icom = 1, 3
            dist = dist + (kvc(icom, 1) - kvc(icom, i))**2
         enddo
         if (dm - dist > eps) then
            i2 = i
            dm = dist
         endif
      enddo
      IF (i2 == 0) CALL juDFT_error(" tetcon1 ", calledby="tetcon")
      dnorm = sqrt(dm)
      ntetra(2, 1) = i2

! FIND POINT NEAREST TO (1+I2)/2 , NOT ON THE LINE
! CONNECTING 1 AND I2

      dm = 10.0**9
      i3 = 0
      do i = 2, nkpt
         if (i == i2) cycle
         dl = 0.0
         dist = 0.0
         do icom = 1, 3
            vect = kvc(icom, i) - 0.5*(kvc(icom, 1) + kvc(icom, i2))
            dist = dist + vect*vect
            dl = dl + vect*(kvc(icom, 1) - kvc(icom, i2))
         enddo
         dl = dl/dnorm
         if (abs(dist - dl*dl) < 0.01*dist) cycle
         if (dm - dist > eps) then
            dm = dist
            i3 = i
         endif
      enddo
      IF (i3 == 0) CALL juDFT_error(" tetcon2 ", calledby="tetcon")
      ntetra(3, 1) = i3

! FIND POINT NEAREST TO (1+I2+I3)/3 , NOT IN THE
! PLANE (1,I2,I3)

      cn(1) = (kvc(2, 1) - kvc(2, i2))*(kvc(3, i2) - kvc(3, i3)) - &
              (kvc(3, 1) - kvc(3, i2))*(kvc(2, i2) - kvc(2, i3))
      cn(2) = (kvc(3, 1) - kvc(3, i2))*(kvc(1, i2) - kvc(1, i3)) - &
              (kvc(1, 1) - kvc(1, i2))*(kvc(3, i2) - kvc(3, i3))
      cn(3) = (kvc(1, 1) - kvc(1, i2))*(kvc(2, i2) - kvc(2, i3)) - &
              (kvc(2, 1) - kvc(2, i2))*(kvc(1, i2) - kvc(1, i3))
      dnorm = 0.0
      do icom = 1, 3
         dnorm = dnorm + cn(icom)**2
      enddo
      dnorm = 1.0/sqrt(dnorm)
      do icom = 1, 3
         cn(icom) = cn(icom)*dnorm
      enddo
      i4 = 0
      dm = 10.0**9
      do i = 2, nkpt
         if ((i == i2) .OR. (i == i3)) cycle
         dist = 0.0
         dl = 0.0
         do icom = 1, 3
            vect = kvc(icom, i) - (kvc(icom, 1) + kvc(icom, i2) + &
                                   kvc(icom, i3))/3.0
            dist = dist + vect**2
            dl = dl + vect*cn(icom)
         enddo
         if (dl*dl < 0.01*dist) cycle
         if (dm - dist > eps) then
            i4 = i
            dm = dist
            vt = dl/(dnorm*6.0)
         endif
      enddo
      IF (i4 == 0) CALL juDFT_error(" tetcon3 ", calledby="tetcon")
      ntetra(4, 1) = i4
      voltet(1) = abs(vt)

! ENTER LOOP FOR CONSTRUCTION OF TETRAHEDRA:

      nt = 1
      it = 0
1000  continue
      it = it + 1

! CHECK THE SIDES OF TETRAHEDRON IT

      do 1100 j = 1, 4

         ! CHECK SIDE OPPOSITE TO CORNER J:

         ic = 0
         do 1200 i = 1, 4
            if (i == j) goto 1200
            ic = ic + 1
            nside(ic) = ntetra(i, it)
            do 1300 icom = 1, 3
               kcorn(icom, ic) = kvc(icom, ntetra(i, it))
1300        END DO
1200     END DO
         nside(4) = ntetra(j, it)
         is1 = min(nside(1), nside(2), nside(3))
         is3 = max(nside(1), nside(2), nside(3))
         is2 = nside(1) + nside(2) + nside(3) - is1 - is3

         ! CHECK IF THERE IS ALREADY A TETRAHEDRON CONNECTED
         ! TO THIS SIDE:

         do nnt = 1, nt
            do 1310 nsid = 1, 4
               if (nnt == it) goto 1310
               ic = 0
               do 1320 i = 1, 4
                  if (i == nsid) goto 1320
                  ic = ic + 1
                  nside2(ic) = ntetra(i, nnt)
1320           END DO
               js1 = min(nside2(1), nside2(2), nside2(3))
               js3 = max(nside2(1), nside2(2), nside2(3))
               js2 = nside2(1) + nside2(2) + nside2(3) - js1 - js3
               if ((is1 == js1) .AND. (is2 == js2) .AND. (is3 == js3)) &
                  goto 1100
1310        END DO
         end do

         ! CONSTRUCT THE OUTWARD NORMAL ON THIS SIDE

         norm(1) = (kcorn(2, 2) - kcorn(2, 1))*(kcorn(3, 3) - kcorn(3, 1)) - &
                   (kcorn(3, 2) - kcorn(3, 1))*(kcorn(2, 3) - kcorn(2, 1))
         norm(2) = (kcorn(3, 2) - kcorn(3, 1))*(kcorn(1, 3) - kcorn(1, 1)) - &
                   (kcorn(1, 2) - kcorn(1, 1))*(kcorn(3, 3) - kcorn(3, 1))
         norm(3) = (kcorn(1, 2) - kcorn(1, 1))*(kcorn(2, 3) - kcorn(2, 1)) - &
                   (kcorn(2, 2) - kcorn(2, 1))*(kcorn(1, 3) - kcorn(1, 1))
         vol = 0.0
         do 1400 i = 1, 3
            vol = vol + norm(i)*(kvc(i, ntetra(j, it)) - kcorn(i, 1))
1400     END DO
         vol = vol/6.0
         if (vol < 0.0) goto 1500
         do 1600 i = 1, 3
            norm(i) = -norm(i)
1600     END DO
1500     continue
         vol = abs(vol)

         ! STORE THE K-POINT ADDRESS IN NKADD ARRAY ACCORDING TO THE
         ! ORDER OF DISTANCE BETWEEN THE K-POINT AND THIS FACE.

         do 1800 nkp = 1, nkpt
            length = 0.0
            do 1850 i = 1, 3
               length = length + (kvc(i, nkp) &
                                  - (kcorn(i, 1) + kcorn(i, 2) + kcorn(i, 3))/3.0)**2
1850        END DO
            ddist(nkp) = length
            nkadd(nkp) = nkp
1800     END DO
         do 1900 nkp = 1, nkpt - 1
            minlen = ddist(nkp)
            ik = nkp
            do 1950 nkq = nkp + 1, nkpt
               if (ddist(nkq) > minlen - eps) goto 1950
               minlen = ddist(nkq)
               ik = nkq
1950        END DO
            ddist(ik) = ddist(nkp)
            k = nkadd(nkp)
            nkadd(nkp) = nkadd(ik)
            nkadd(ik) = k
1900     END DO

         ! CONSTRUCT A TETRAHEDRON WHICH CONNECT THIS FACE TO
         ! A K-POINT.

         ik = 0
         do 2000 nkp = 1, nkpt
            do 2010 i = 1, 3
               if (nkadd(nkp) == nside(i)) goto 2000
2010        END DO
            do 2050 i = 1, 3
               kcorn(i, 4) = kvc(i, nkadd(nkp))
2050        END DO
            nside(4) = nkadd(nkp)
            vt = 0.0
            do 2100 i = 1, 3
               vt = vt + norm(i)*(kcorn(i, 4) - kcorn(i, 1))
2100        END DO
            vt = vt/6.0

            ! REJECT POINT NKP IF IT IS ON THE WRONG SIDE .

            if (vt < eps*vol) goto 2000

            ! CHECK IF THIS TETRAHEDRON INTERSECTS AN EXISTING ONE

            do 3100 nk2 = 1, nt
               if (nk2 == it) goto 3100

               ! FIRST CHECK IF TETRAHEDRON NK2 IS ON THE DANGEROUS SIDE

               do 3150 i = 1, 4
                  sum = 0.0
                  do 3160 jj = 1, 3
                     sum = sum + norm(jj)*(kvc(jj, ntetra(i, nk2)) - kcorn(jj, 1))
3160              END DO
                  if (sum > eps) goto 3170
3150           END DO
               goto 3100
3170           continue

               ! TETRAHEDRON NK2 IS POTENTIALLY SUSPECT

               l = 0
               do i = 1, 4
                  do 3200 jj = 1, 4
                     if (ntetra(i, nk2) == nside(jj)) goto 3200
                     sum = 0.0
                     do icom = 1, 3
                        cn(icom) = kvc(icom, ntetra(i, nk2)) - kcorn(icom, jj)
                        sum = sum + cn(icom)*cn(icom)
                     end do
                     do icom = 1, 3
                        cn(icom) = cn(icom)/sqrt(sum)
                     end do
                     do k = 1, l
                        sum = 0.0
                        do icom = 1, 3
                           sum = sum + cn(icom)*d(icom, k)
                        end do
                        if ((1.0 - sum) < eps) goto 3200

                        ! EXCLUDE THE VECTOR WHICH HAS THE SAME DIRECTION AS AN EXISTING ONE.

                     end do
                     l = l + 1
                     do icom = 1, 3
                        d(icom, l) = cn(icom)
                     end do
3200              END DO
               end do
               IF (l < 4) CALL juDFT_error(" tetcon9 ", calledby="tetcon")

               ! HERE, WE HAVE A SET OF D-VECTORS WHICH CONNECT THE CORNER POINTS
               ! OF ONE TETRAHEDRON NK2 WITH THOSE OF THE CURRENT TETRAHEDRON TO
               ! BE CHECKED.
               ! FIND A PLANE OF WHICH ALL D-VECTORS EXIST IN ONE SIDE.

               do 3510 i = 1, l - 1
                  do 3500 jj = i + 1, l
                     dnm(1) = d(2, i)*d(3, jj) - d(3, i)*d(2, jj)
                     dnm(2) = d(3, i)*d(1, jj) - d(1, i)*d(3, jj)
                     dnm(3) = d(1, i)*d(2, jj) - d(2, i)*d(1, jj)

                     ! DNM IS THE NORMAL VECTOR TO THE PLANE GIVEN BY
                     ! I-TH AND JJ-TH D-VECTORS.

                     if ((abs(dnm(1)) < eps) .AND. (abs(dnm(2)) < eps) &
                         .AND. (abs(dnm(3)) < eps)) goto 3500
                     do 3600 k = 1, l
                        if ((k == i) .OR. (k == jj)) goto 3600
                        sum = 0.0
                        do icom = 1, 3
                           sum = sum + d(icom, k)*dnm(icom)
                        end do
                        if (abs(sum) < eps) goto 3600
                        do 3700 kk = 1, l
                           if ((kk == i) .OR. (kk == jj) .OR. (kk == k)) goto 3700
                           sss = 0.0
                           do icom = 1, 3
                              sss = sss + d(icom, kk)*dnm(icom)
                           end do
                           if (abs(sss) < eps) goto 3700
                           if (sss*sum < 0.0) goto 3500

                           ! IF K-TH AND KK-TH D-VECTORS EXIST IN OPPOSITE SIDE WITH
                           ! RESPECT TO THE PLANE GIVEN BY I-TH AND JJ-TH D-VECTORS,
                           ! WE WILL ATTEMPT THE NEXT PLANE. (GOTO 3500)

3700                    END DO
                        goto 3100

                        ! WE SUCCEED WHEN THE PLANE SATISFIES THE CONDITION,
                        ! I.E., THE CURRENT TETRAHEDRON DOES NOT INTERSECTS NK2-TH ONE:
                        ! WE TRY TO CHECK THE NEXT TETRAHEDRON. (GOTO 3100)

3600                 END DO
3500              END DO
3510           END DO
               goto 2000

               ! HERE, A TETRAHEDRON MADE OF J-TH FACE AND NKADD(NKP) K-POINT
               ! INTERSECTS AT LEAST ONE OF THE EXISTING TETRAHEDRONS:
               ! WE TRY THE NEXT K-POINT. (GOTO 2000)

3100        END DO
            ik = nkadd(nkp)
            volnew = vt
            goto 2500

            ! HERE, WE GET THE NEW TETRAHEDRON WHICH IS MADE OF J-TH FACE
            ! AND IK-TH K-POINT.

2000     END DO
         goto 1100

         ! HERE, WE COULD NOT FIND ANY NEW TETRAHEDRON USING J-TH FACE.
         ! WE TRY THE NEXT FACE. (GOTO 1100)

2500     continue
         nt = nt + 1
         do 2300 i = 1, 4
            ntetra(i, nt) = ntetra(i, it)
2300     END DO
         ntetra(j, nt) = ik
         voltet(nt) = abs(volnew)

1100  END DO
      if (it < nt) goto 1000
      IF (nt > ndiv3) CALL juDFT_error(" nt>ndiv3", calledby="tetcon")

! DETERMINE THE CHARACTERISTICS OF THIS DIVISION INTO
! TETRAHEDRA .

      vav = 0.0
      vsq = 0.0
      sav = 0.0
      ssq = 0.0
      vmin = 10.0**9
      vmax = -vmin
      smin = vmin
      smax = vmax
      do it = 1, nt
         vt = voltet(it)
         vav = vav + vt
         vsq = vsq + vt**2
         if (vt > vmax) vmax = vt
         if (vt < vmin) vmin = vt
         do i = 1, 3
            do j = i + 1, 4
               dist = 0.0
               do icom = 1, 3
                  dist = dist + (kvc(icom, ntetra(i, it)) - &
                                 kvc(icom, ntetra(j, it)))**2
               end do
               dist = sqrt(dist)
               sav = sav + dist
               ssq = ssq + dist*dist
               if (dist > smax) smax = dist
               if (dist < smin) smin = dist
            end do
         end do
      end do
      vav = vav/nt
      sav = sav/(6*nt)
      vsq = vsq/nt
      ssq = ssq/(6*nt)
      vsq = sqrt(vsq - vav*vav)
      ssq = sqrt(ssq - sav*sav)
      write (oUnit, 5000) nt, vav, vsq, vmin, vmax, sav, ssq, smin, smax
5000  format(/, '   division into tetrahedra  ', /, &
              '  there are      ', i5, '  tetrahedra ', /, &
              '  volume         ', f15.10, '  +/-  ', f10.5, 3x, 2f10.5, /, &
              '  side           ', f15.10, '  +/-  ', f10.5, 3x, 2f10.5,/)
!     write(oUnit,5000) nt,vav,vsq,vmin,vmax,sav,ssq,smin,smax
      write (oUnit, 5100) ((ntetra(j, i), j=1, 4), i=1, nt)
5100  format(4(4x, 4i4))
!     write(oUnit,5100) ((ntetra(j,i),j=1,4),i=1,nt)
! CHECK IF WE HAVE THE CORRECT TOTAL VOLUME
      vt = omega*vav*nt*nsym/(2*pi_const)**3 - 1.0
      write (oUnit, 5200) vt
!     write(oUnit,5200) vt
      do 5300 i = 1, nt
         !     write(oUnit,'('' tetrahedra # '',i5,'' is '',d12.4)') i,voltet(i)
5300  END DO
5200  format(/, '  voltetsum/volBZ - 1  ', d12.4)

!     The following statement used to have a stop in it.
!     If the word TETCON5 appears you have failed the < 1.0D-5 test.

      if (abs(vt) > eps1) THEN
         write (oUnit, '(''  tetcon5  '')')
      endif
      RETURN
   END SUBROUTINE tetcon
END MODULE m_tetcon
