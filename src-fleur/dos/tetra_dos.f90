!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_tetra_dos
   !----------------------------------------------------------------------
   !
   ! This subroutine evaluates the density of states (g) by the linear
   ! tetrahedron method on a energy grid (e) of 'ned' points.
   !
   ! eig()          ... eigenvalues
   ! qal()         ... partial charges
   ! ntet          ... number of tetrahedrons
   ! ntetra(1-4,ntet)... index of k-points forming tetrahedron nt
   ! voltet(ntet)    ... volume of tetrahedron nt
   !
   !                                                      gb 2000
   !----------------------------------------------------------------------
   USE m_types_kpts

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE tetra_dos(jspins,kpts,eGrid,neig,eig,qal,g,energyShift)

      INTEGER,       INTENT(IN)    :: jspins
      TYPE(t_kpts),  INTENT(IN)    :: kpts
      INTEGER,       INTENT(IN)    :: neig(:,:)     !(nkpt,jspins)
      REAL,          INTENT(IN)    :: eGrid(:)   !(ned)
      REAL,          INTENT(IN)    :: qal(:,:,:)  !(neigd,nkpt,jspins)
      REAL,          INTENT(IN)    :: eig(:,:,:)  !(neigd,nkpt,jspins)
      REAL,          INTENT(OUT)   :: g(:,:)      !(ned,jspins)
      REAL, OPTIONAL, INTENT(IN)   :: energyShift

      INTEGER :: i,j,iBand,ikpt,ie,idim,itet,icorn,jcorn,ispin
      REAL    :: ener,w,sfac, shift, term
      REAL    :: weight(4),eval(4),ecmax(SIZE(eig,1),jspins)
      REAL, ALLOCATABLE :: wpar(:,:,:)
      REAL, ALLOCATABLE :: eig_nondeg(:,:,:)

      shift = 0.0
      IF(PRESENT(energyShift)) shift = energyShift

      eig_nondeg = eig - shift

      ALLOCATE(wpar,mold=qal)
      wpar = 0.0

      sfac = 2.0/jspins

      DO ispin = 1, SIZE(qal,3)
         DO iBand = 1,SIZE(eig,1)
            ecmax(iBand,ispin) = -1.0e25
            DO ikpt = 1,kpts%nkpt
               IF(eig_nondeg(iBand,ikpt,ispin).GT.ecmax(iBand,ispin)) ecmax(iBand,ispin) = eig_nondeg(iBand,ikpt,ispin)
            ENDDO
         ENDDO
      ENDDO
      !
      !  check for energy degeneracies in tetrahedrons
      !
      DO ispin = 1, SIZE(qal,3)
         DO itet = 1,kpts%ntet
            DO iBand = 1,SIZE(eig_nondeg,1)
               DO i = 1,3
                  icorn = kpts%ntetra(i,itet)
                  DO j = i+1,4
                     jcorn = kpts%ntetra(j,itet)
                     IF (abs(eig_nondeg(iBand,icorn,ispin)-eig_nondeg(iBand,jcorn,ispin)).LT.1.0e-7) THEN
                        eig_nondeg(iBand,icorn,ispin) = eig_nondeg(iBand,icorn,ispin) + i*1.0e-7*itet
                        eig_nondeg(iBand,jcorn,ispin) = eig_nondeg(iBand,jcorn,ispin) - i*1.0e-7*itet
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      !WRITE (*,*) 'in tetra_dos'  ! do not remove  this statement

      !
      ! calculate partial weights
      !
      DO ispin = 1, SIZE(qal,3)
         DO ikpt=1,kpts%nkpt
            DO iBand = 1,neig(ikpt,ispin)
               DO itet = 1,kpts%ntet
                  IF (ALL(kpts%ntetra(:,itet).ne.ikpt)) CYCLE

                  eval = eig_nondeg(iBand,kpts%ntetra(:,itet),ispin)

                  IF(ANY(eval.GE.9999.9)) CYCLE

                  DO i=1,4
                     icorn = kpts%ntetra(i,itet)

                     weight(i)=1.0
                     DO j=1,4
                        IF (i.NE.j) weight(i)=weight(i)*(eval(j)-eval(i))
                     ENDDO
                     weight(i)=6.0*kpts%voltet(itet)/(weight(i)*kpts%ntet)

                     wpar(iBand,icorn,ispin) =  wpar(iBand,icorn,ispin) &
                                             + sfac*0.25*weight(i)*qal(iBand,ikpt,ispin)
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      !---------------------------------------------------
      ! evaluate d.o.s.
      !---------------------------------------------------
      !
      g = 0.0

      DO ispin = 1, SIZE(qal,3)
         DO ikpt = 1,kpts%nkpt
            DO iBand = 1,neig(ikpt,ispin)

               IF(MINVAL(eGrid).GT.ecmax(iband,ispin)) cycle

               ener = eig_nondeg(iBand,ikpt,ispin)
               w  = 0.5*wpar(iBand,ikpt,ispin)
               DO ie = 1,SIZE(eGrid)
                  term = eGrid(ie) - ener
                  IF(eGrid(ie).GT.ecmax(iBand,ispin)) cycle
                  IF(term.LT.0.0e0) cycle
                  g(ie,ispin) = g(ie,ispin) + w * term**2
               ENDDO

            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE tetra_dos
END MODULE m_tetra_dos
