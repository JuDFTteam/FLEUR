MODULE m_tetra

   !This module contains subroutines used in the DFT+HUBBARD1 method to calculate the 
   !imaginary part of the on-site green's function using the linear tetrahedron method 
   !used in an old version of KKR

CONTAINS

SUBROUTINE int_tetra(nkpt,ntetra,itetra,voltet,neigd,nevk,qal,eig,ne,g,top,bot)

   IMPLICIT NONE

   INTEGER,                INTENT(IN)     :: ntetra, nkpt,ne
   REAL,                   INTENT(IN)     :: voltet(ntetra)
   INTEGER,                INTENT(IN)     :: itetra(4,6*nkpt)

   INTEGER,                INTENT(IN)     :: neigd
   INTEGER,                INTENT(IN)     :: nevk(nkpt)
   REAL,                   INTENT(IN)     :: qal(neigd,nkpt)
   REAL,                   INTENT(IN)     :: eig(neigd,nkpt)
   REAL   ,                INTENT(INOUT)  :: g(ne)
   REAL,                   INTENT(IN)     :: top
   REAL,                   INTENT(IN)     :: bot

   INTEGER ntet, i, j, k, l
   REAL tmp, vol
   LOGICAL valid
   REAL mean, e(4)

   DO ntet = 1, ntetra

      valid = .true.
      !calculate the mean value (valid?) of qal inside one tetrahedron and order the eigenenergies
      vol = voltet(ntet)/ntetra

      DO i = 1, neigd

         mean = 0.0
         DO j = 1, 4
            mean = mean + 0.25 * qal(i,itetra(j,ntet))
         ENDDO

         DO j = 1, 4
            IF (i.GT.nevk(itetra(j,ntet))) THEN !Check wether the eigenvalue with band index i 
                                                !exist at every corner 
               valid = .false.
            ELSE
               e(j) = eig(i,itetra(j,ntet))
            ENDIF
         ENDDO

         IF(valid) THEN
            !Order from lowest to highest
            DO l = 1, 3
               DO k = l+1, 4
                  IF (e(l).GT.e(k)) THEN
                     tmp = e(l)
                     e(l) = e(k)
                     e(k) = tmp
                  ENDIF
               ENDDO
            ENDDO

            !check for energy degeneracies (to avoid nan in weight factors)
            DO l = 1, 3
               DO k = l+1, 4
                  IF(abs(e(l)-e(k)).LT.1.0E-9) THEN
                     e(l) = e(l) - i*1.0E-9
                     e(k) = e(k) + i*1.0E-9
                  ENDIF
               ENDDO
            ENDDO

            CALL contr_singletetra(g,mean,e,vol,ne,top,bot)
         ENDIF
      ENDDO

   ENDDO


END SUBROUTINE int_tetra

SUBROUTINE contr_singletetra(g,mean,ev,vol,ne,top,bot)

   USE m_constants

   IMPLICIT NONE

   INTEGER,                INTENT(IN)     :: ne
   REAL,                   INTENT(INOUT)  :: ev(4)
   REAL,                   INTENT(IN)     :: mean
   REAL,                   INTENT(IN)     :: vol, top, bot
   REAL,                   INTENT(INOUT)  :: g(ne)
   
   INTEGER i , j, n
   REAL emin, emax, e
   REAL weight 
   REAL a,b,c
   INTEGER ind(3)
   REAL del, deln

   REAL aw, bw
   INTEGER nstart, nend


   del = (top-bot)/REAL(ne-1)
   ind(1:3) = (/1,2,4/)
   !check if the tetrahedron is insíde the desired energy range to avoid calculating 
   !tetrahedrons completely outside of that range
   IF((ev(1).LT.top).AND.(ev(4).GT.bot)) THEN
      DO i = 1, 3
         emin = ev(i)
         emax = ev(i+1)

         !determine the last index with a lower energy than emin/max
         nstart = INT((emin-bot)/del)+1
         nend = INT((emax-bot)/del)+1

         !distance between e(nstart+1) and emin 
         aw = del - (emin - ((nstart-1)*del+bot))
         !distance between e(nend) and emax
         bw = emax - ((nend-1)*del+bot)

         weight = 3.0 * vol
         DO j = 1, 4
            IF (ind(i).NE.j) weight = weight/ABS(ev(ind(i))-ev(j))
         ENDDO

         !calculate the right prefactors
         IF(i.EQ.2) THEN
            a = a - emin**2 * weight
            b = b + emin * weight
            c = c - 1./3. * weight
         ELSE
            IF(i.EQ.3) THEN
               a = emax**2 * weight
               b = (-1.0) * emax * weight
               c = 1./3. * weight
            ELSE 
               a = emin**2 * weight
               b = (-1.0) * emin * weight
               c = 1./3. * weight
            ENDIF
         ENDIF

         !If emin and emax are inside the same energy step we use a different formula (from where?)
         IF(nstart.EQ.nend) THEN
            deln = (emax-emin)/del * (a + b*(emax+emin) + c *((emax+emin)**2-emax*emin))
            IF((nstart.LE.ne).AND.(nstart.GE.1)) g(nstart) = g(nstart) + deln * mean
         ELSE
            !Part of the integral between emin and the first grid point
            deln = aw/del * (a + b*(2*emin+aw) + c*((2*emin+aw)**2-emin*(emin+aw)))
            IF((nstart.LE.ne).AND.(nstart.GE.1)) g(nstart) = g(nstart) + deln * mean
            e = emin + aw
            !loop over grid points between emin and emax
            DO n = nstart + 1, nend-1 
               deln = (a + b*(2*e+del) + c*((2*e+del)**2-e*(e+del)))
               IF((n.LE.ne).AND.(n.GE.1)) g(n) = g(n) + deln * mean
               e = e + del
            ENDDO
            !Part of the integral between the last grid point before emax and emax
            deln = bw/del * (a + b*(2*e+bw) + c*((2*e+bw)**2-e*(e+bw)))
            IF((nend.LE.ne).AND.(nend.GE.1)) g(nend) = g(nend) + deln * mean
         ENDIF
      ENDDO
   ENDIF  



END SUBROUTINE contr_singletetra





END MODULE m_tetra