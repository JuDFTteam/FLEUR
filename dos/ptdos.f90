MODULE m_ptdos
   !-------------------------------------------------------------------------
   ! Density of states calculated by linear triangular method
   !-------------------------------------------------------------------------
   IMPLICIT NONE

   CONTAINS

   SUBROUTINE ptdos(jspins,ne,ndos,ntb,ntria,as,atr,&
                       itria,nkpt,ev,qal,e,g)

      INTEGER, INTENT (IN) :: ne,ntria,jspins,ndos,ntb,nkpt
      INTEGER, INTENT (IN) :: itria(:,:)
      REAL,    INTENT (IN) :: as
      REAL,    INTENT (IN) :: atr(:),qal(:,:,:)!(ndos,ntb,nkpt)
      REAL,    INTENT (IN) :: e(:),ev(:,:)!(ntb,nkpt)
      REAL,    INTENT (OUT):: g(:,:)      !(ne,ndos)

      INTEGER :: i, j, nl, nb, n, nt(3), nc(4)
      REAL    :: f, fa, ec(4)

      !
      !     calculate partial densities of states
      !
      f = 2*(3-jspins)/as

      g = 0.

      DO n = 1 , ntria
         fa = f*atr(n)
         nt(:) = itria(:,n)
         DO nb = 1 , ntb
            ec(1:3) = ev(nb,nt(:))
            nc(1:3) = nt(:)
            DO i = 1, 2
               DO j = i+1, 3
                  IF ( ec(i).GT.ec(j) ) THEN
                     ec(4) = ec(i) ; ec(i) = ec(j) ; ec(j) = ec(4)
                     nc(4) = nc(i) ; nc(i) = nc(j) ; nc(j) = nc(4)
                  ENDIF
               ENDDO
            ENDDO

            DO nl = 1 , ndos
               DO i = 1 , ne
                  g(i,nl) = g(i,nl) + fa* dostet( e(i),ec(1),ec(2),ec(3),&
                              qal(nl,nb,nc(1)),qal(nl,nb,nc(2)),qal(nl,nb,nc(3)) )
               ENDDO
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE ptdos

   REAL FUNCTION dostet(e,e1,e2,e3,q1,q2,q3)
      !
      !     partial density of states for one tetrahedron
      !     note that e1.le.e2.le.e3 is assumed
      !

      REAL, INTENT(IN) :: e , e1 , e2 , e3 , q1 , q2 , q3

      REAL :: e21 , e31 , e32 , ee
      REAL, PARAMETER :: tol = 1.e-6

      dostet = 0.
      IF ( e.LT.e1 ) RETURN
      IF ( e.LE.e2 ) THEN
         !
         !     case 1: e between e1 and e2
         !
         ee = e - e1
         e21 = e2 - e1
         IF ( e21.LT.tol ) RETURN
         e31 = e3 - e1
         IF ( e31.LT.tol ) RETURN
         dostet = ee/(e21*e31)*(q1+0.5*(ee/e21*(q2-q1)+ee/e31*(q3-q1)))
         RETURN
      ELSE
         IF ( e.GT.e3 ) RETURN
         !
         !     case 2: e between e2 and e3
         !
         e31 = e3 - e1
         IF ( e31.LT.tol ) RETURN
         e32 = e3 - e2
         IF ( e32.LT.tol ) RETURN
         dostet = (e3-e)/(e31*e32)&
                  *0.5*(q1+q2+(e-e1)/e31*(q3-q1)+(e-e2)/e32*(q3-q2))
         RETURN
      ENDIF
   END FUNCTION dostet

END MODULE m_ptdos
