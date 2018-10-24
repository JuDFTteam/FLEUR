      MODULE m_triang
      use m_juDFT
!-------------------------------------------------------------------
c     find a triangular decomposition of the irreducible wedge of
c     the first brillouin zone for a given k-mesh. k-points at all
c     vertices are assumed.
c     erich wimmer     july 1981
!-------------------------------------------------------------------

      IMPLICIT NONE

      CONTAINS
      SUBROUTINE triang(
     >                  v,nkpt,
     <                  it,ntria,at,att,l_f_t)

c     Arguments
      INTEGER, INTENT(IN)  :: nkpt
      REAL,    INTENT(IN)  :: v(3,nkpt)
      REAL,    INTENT(OUT) :: att , at(2*nkpt)
      INTEGER, INTENT(OUT) :: ntria , it(3,2*nkpt)
      LOGICAL, INTENT(INOUT) :: l_f_t
c     locals
      REAL    :: a, a1, d, dm, s0, s1, x1, x2, x3, y1, y2, y3
      INTEGER :: j , j1 , j2 , jc , jj , k , k1 , k2 , k3 , kk , l1 , 
     &           l2 , l3 , n , nt0 , nkpts
      LOGICAL :: new
c     constants
      REAL , PARAMETER :: zero = 0.0 , big = 1.e8 , tol = 1.e-5

      ntria = 0
      IF ( nkpt.LT.3 ) RETURN
c
c l_f_t=.true. means that we call from fertri and on output gives 'film'
c
      IF (l_f_t) THEN ! check for bulk
         DO j= 2,nkpt
            IF (abs(v(3,j)-v(3,1)).GT.1.0e-12) THEN
              l_f_t = .false.
              nkpts = j - 1
              IF ((mod(nkpt,nkpts).NE.0).OR.(j.LT.3)) THEN
                WRITE (6,*) 'tria=T & film=F requires k-point planes'
                WRITE (6,*) 'with equally distributed k-points !'
c                 CALL juDFT_error("not a k-point set for tria=T & film=F",calledby="triang")
              ENDIF
!              RETURN
              GOTO 10
            ENDIF
         ENDDO
         l_f_t = .true.
         nkpts = nkpt
 10      CONTINUE
      ELSE
c
c l_f_t=.false. means that we call from evaldos
c
         nkpts = nkpt
         DO j= 2,nkpt
            IF (abs(v(3,j)-v(3,1)).GT.1.0e-12) THEN
              nkpts = j - 1
              GOTO 11
           ENDIF
        ENDDO
 11     CONTINUE
        IF ((mod(nkpt,nkpts).NE.0).OR.(j.LT.3)) THEN
          l_f_t = .false.
          RETURN
        ELSE
          l_f_t = .true.
        ENDIF
      ENDIF
c
c     create first triangle
c
      att = zero
      k1 = 1
      dm = big
      DO k = 2 , nkpts
         d = vd2(v(1,K1),v(2,k1),v(1,k),v(2,k))
         IF ( d.LT.dm ) THEN
            dm = d
            k2 = k
         ENDIF
      ENDDO
c
      dm = big
      new = .false.
      DO k = 1 , nkpts
         IF ( k.NE.k1 .AND. k.NE.k2 ) THEN
            d = VD2(V(1,k1),V(2,k1),V(1,k),V(2,k))
     &          + VD2(V(1,k2),V(2,k2),V(1,k),V(2,k))
            a = AREA(V(1,k1),V(2,k1),V(1,k2),V(2,k2),V(1,k),V(2,k))
            IF ( ABS(a).GE.TOL ) THEN
               IF ( d.LT.dm ) THEN
                  new = .TRUE.
                  dm = d
                  a1 = a
                  k3 = k
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF ( new ) THEN
         Ntria = 1
         It(1,Ntria) = k1
         IF ( k3.LE.k2 ) THEN
            kk = k2
            k2 = k3
            k3 = kk
         ENDIF
         It(2,Ntria) = k2
         It(3,Ntria) = k3
         At(Ntria) = ABS(a1)/2
         Att = Att + At(Ntria)
c
c     create, if possible, a new triangle from each side of the mother
c     triangle nt0 with minimal sum of sides
c
         nt0 = 0
         DO WHILE ( .TRUE. )
            nt0 = nt0 + 1
            IF ( nt0.GT.Ntria ) RETURN
            DO l1 = 1 , 3
               l2 = MOD(l1,3) + 1
               l3 = MOD(l2,3) + 1
               k1 = It(l1,nt0)
               k2 = It(l2,nt0)
               IF ( k2.LE.k1 ) THEN
                  kk = k1
                  k1 = k2
                  k2 = kk
               ENDIF
c--->>     check if side k1-k2 belongs already to another triangle
               new = .TRUE.
               DO n = 1 , Ntria
                  IF ( n.NE.nt0 ) THEN
                     IF ( (k1.EQ.It(1,n) .AND. k2.EQ.It(2,n)) .OR. 
     &                    (k1.EQ.It(1,n) .AND. k2.EQ.It(3,n)) .OR. 
     &                    (k1.EQ.It(2,n) .AND. k2.EQ.It(3,n)) )
     &                    new = .FALSE.
                  ENDIF
               ENDDO
               IF ( new ) THEN
                  k3 = It(l3,nt0)
                  a = AREA(V(1,k1),V(2,k1),V(1,k2),V(2,k2),V(1,k3),
     &                V(2,k3))
                  s0 = SIGN(1.,a)
c--->>     a new triangle sharing the side k1-k2 with the mother
c--->>     triangle nt0 has to ly on the other side, i.e. the cross
c--->>     products (k2-k1)x(k3(old)-k1) and (k2-k1)x(k3(new)-k1)
c--->>     have to have opposite signs
                  dm = BIG
                  new = .FALSE.
                  DO k = 1 , Nkpts
                     IF ( k.NE.k1 .AND. k.NE.k2 ) THEN
c--->>     check if a new side, (k1,k) or (k2,k), belongs
c--->>     already to an older triangle
                        j1 = k1
                        j2 = k
                        DO j = 1 , 2
                           IF ( j2.LE.j1 ) THEN
                              jj = j1
                              j1 = j2
                              j2 = jj
                           ENDIF
                           jc = 0
                           DO n = 1 , Ntria
                              IF ( j1.EQ.It(1,n) .AND. j2.EQ.It(2,n)
     &                             .OR. j1.EQ.It(1,n) .AND. 
     &                             j2.EQ.It(3,n) .OR. j1.EQ.It(2,n)
     &                             .AND. j2.EQ.It(3,n) ) jc = jc + 1
                           ENDDO
                           IF ( jc.EQ.2 ) GOTO 5
                           j1 = k2
                           j2 = k
                        ENDDO
                        a = AREA(V(1,k1),V(2,k1),V(1,k2),V(2,k2),V(1,k),
     &                      V(2,k))
                        IF ( ABS(a).GE.TOL ) THEN
                           s1 = SIGN(1.,a)
                           IF ( ABS(s1-s0).GE.TOL ) THEN
                              d = VD2(V(1,k1),V(2,k1),V(1,k),V(2,k))
     &                            + VD2(V(1,k2),V(2,k2),V(1,k),V(2,k))
                              IF ( d.LT.dm ) THEN
                                 new = .TRUE.
                                 dm = d
                                 a1 = a
                                 k3 = k
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
 5                ENDDO
                  IF ( new ) THEN
                     Ntria = Ntria + 1
                     IF ( Ntria>2*nkpt )  CALL juDFT_error("ntriad"
     +                    ,calledby ="triang")
c
                     IF ( k2.LE.k1 ) THEN
                        kk = k1
                        k1 = k2
                        k2 = kk
                     ENDIF
                     IF ( k3.LE.k1 ) THEN
                        kk = k1
                        k1 = k3
                        k3 = kk
                     ENDIF
                     IF ( k3.LE.k2 ) THEN
                        kk = k2
                        k2 = k3
                        k3 = kk
                     ENDIF
                     It(1,Ntria) = k1
                     It(2,Ntria) = k2
                     It(3,Ntria) = k3
                     At(Ntria) = ABS(a1)/2
                     Att = Att + At(Ntria)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE
c     write(16,1000) ((v(i,j),i=1,2),j=1,nt)
!          CALL juDFT_error("triang",calledby="triang")
      ENDIF
c
99001 FORMAT (' $$$ error in triang: collinear k-points'/(5x,2F12.6))
c
      END SUBROUTINE triang
!-----------------------------------------------------
      REAL FUNCTION vd2(x1,y1,x2,y2)         ! distance between (x1,y1) and (x2,y2)
         REAL x1,y1,x2,y2
         vd2 = (x2-x1)**2 + (y2-y1)**2
      END FUNCTION vd2
      REAL FUNCTION area(x1,y1,x2,y2,x3,y3)  ! area of triangle (x1,y1),(x2,y2),(x3,y3)
         REAL x1,y1,x2,y2,x3,y3
         area = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
      END FUNCTION area
!-----------------------------------------------------
      END MODULE m_triang
