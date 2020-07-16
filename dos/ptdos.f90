!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_ptdos
   !-------------------------------------------------------------------------
   ! Density of states calculated by linear triangular method
   !-------------------------------------------------------------------------
   USE m_types
   USE m_tetsrt

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE ptdos(jspins,ne,ndos,nbands,kpts,ev,qal,eGrid,g)

      INTEGER,       INTENT(IN)  :: ne,jspins,ndos,nbands
      TYPE(t_kpts),  INTENT(IN)  :: kpts
      REAL,          INTENT(IN)  :: qal(:,:,:)  !(ndos,nbands,nkpt)
      REAL,          INTENT(IN)  :: eGrid(:)    !(ne)
      REAL,          INTENT(IN)  :: ev(:,:)     !(nbands,nkpt)
      REAL,          INTENT(OUT) :: g(:,:)      !(ne,ndos)

      INTEGER :: idos,iBand,itria,iGrid
      INTEGER :: ind(3),k(3)
      REAL    :: sfac,fa
      REAL    :: ei(3)

      !Spin-degeneracy factor
      sfac = 2.0*(3.0-jspins)

      g = 0.0

      DO itria = 1 , kpts%ntet
         fa = sfac*kpts%voltet(itria)/kpts%ntet
         k = kpts%ntetra(:,itria)
         DO iBand = 1 , nbands
            !---------------------------
            !eigenvalues at the corners
            !of the current triangle
            !---------------------------
            ei = ev(iBand,k)

            !sort in ascending order
            ind = tetsrt(ei)

            !
            !calculate partial densities of states
            !

            DO idos = 1 , ndos
               DO iGrid = 1 , ne
                  g(iGrid,idos) = g(iGrid,idos) + fa &
                                  * dostet( eGrid(iGrid),ei(ind),qal(idos,iBand,k(ind)) )
               ENDDO
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE ptdos

   REAL FUNCTION dostet(e,ei,q)
      !
      !     partial density of states for one tetrahedron
      !     note that e1.le.e2.le.e3 is assumed
      !
      !     ei : eigenvalues at the corners of the triangle (sorted)
      !     q : partial charges at the corners

      REAL, INTENT(IN) :: e      !Energy point where the DOS is calculated
      REAL, INTENT(IN) :: ei(:)  !(3)
      REAL, INTENT(IN) :: q(:)   !(3)

      REAL :: e21 , e31 , e32 , ee
      REAL, PARAMETER :: tol = 1e-6

      dostet = 0.0
      IF ( e.LT.ei(1) ) RETURN
      IF ( e.LE.ei(2) ) THEN
         !
         !     case 1: e between e1 and e2
         !
         ee = e - ei(1)
         e21 = ei(2) - ei(1)
         IF ( e21.LT.tol ) RETURN
         e31 = ei(3) - ei(1)
         IF ( e31.LT.tol ) RETURN
         dostet = ee/(e21*e31)*(q(1) &
                                +0.5*(ee/e21*(q(2)-q(1)) &
                                + ee/e31*(q(3)-q(1))))
      ELSE
         IF ( e.GT.ei(3) ) RETURN
         !
         !     case 2: e between e2 and e3
         !
         e31 = ei(3) - ei(1)
         IF ( e31.LT.tol ) RETURN
         e32 = ei(3) - ei(2)
         IF ( e32.LT.tol ) RETURN
         dostet = (ei(3)-e)/(e31*e32)*0.5*(q(1)+q(2) &
                                        +(e-ei(1))/e31*(q(3)-q(1)) &
                                        +(e-ei(2))/e32*(q(3)-q(2)))
         RETURN
      ENDIF
   END FUNCTION dostet

END MODULE m_ptdos
