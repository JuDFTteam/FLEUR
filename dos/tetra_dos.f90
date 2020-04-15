!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_tetrados
   !----------------------------------------------------------------------
   !
   ! This subroutine evaluates the density of states (g) by the linear
   ! tetrahedron method on a energy grid (e) of 'ned' points.
   !
   ! ev()          ... eigenvalues
   ! qal()         ... partial charges
   ! ntet          ... number of tetrahedrons
   ! ntetra(1-4,ntet)... index of k-points forming tetrahedron nt
   ! voltet(ntet)    ... volume of tetrahedron nt
   !
   !                                                      gb 2000
   !----------------------------------------------------------------------
   USE m_types

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE tetra_dos(qdim,neigd,ned,kpts,efermi,energy,nevk,ev,qal,g)

      INTEGER,       INTENT(IN)    :: neigd,ned,qdim
      REAL,          INTENT(IN)    :: efermi
      TYPE(t_kpts),  INTENT(IN)    :: kpts
      INTEGER,       INTENT(IN)    :: nevk(:)     !(nkpt)
      REAL,          INTENT(IN)    :: energy(:)   !(ned)
      REAL,          INTENT(IN)    :: qal(:,:,:)  !(lmax*ntype+3,neigd,nkpt)
      REAL,          INTENT(INOUT) :: ev(:,:)     !(neigd,nkpt)
      REAL,          INTENT(OUT)   :: g(:,:)      !(ned,lmax*ntype+3)

      INTEGER :: i,j,iBand,ikpt,ie,idim,itet,icorn,jcorn
      REAL    :: ener,w
      REAL    :: weight(4),eval(4),ecmax(neigd),term(ned)
      REAL    :: wpar(qdim,neigd,kpts%nkpt)

      DO ikpt = 1,kpts%nkpt
         ev(nevk(ikpt)+1:neigd,ikpt) = 1.0e10
      ENDDO

      wpar = 0.0

      DO iBand = 1,neigd
         ecmax(iBand) = -1.0e25
         DO ikpt = 1,kpts%nkpt
            IF(ev(iBand,ikpt).GT.ecmax(iBand)) ecmax(iBand) = ev(iBand,ikpt)
         ENDDO
      ENDDO
      !
      !  check for energy degeneracies in tetrahedrons
      !
      DO itet = 1,kpts%ntet
         DO iBand = 1,neigd
            DO i = 1,3
               icorn = kpts%ntetra(i,itet)
               DO j = i+1,4
                  jcorn = kpts%ntetra(j,itet)
                  IF (abs(ev(iBand,icorn)-ev(iBand,jcorn)).LT.1.0e-7) THEN
                     ev(iBand,icorn) = ev(iBand,icorn) + i*1.0e-7*itet
                     ev(iBand,jcorn) = ev(iBand,jcorn) - i*1.0e-7*itet
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      WRITE (*,*) 'in tetra_dos'  ! do not remove  this statement

      !
      ! calculate partial weights
      !
      DO ikpt=1,kpts%nkpt
         DO iBand = 1,nevk(ikpt)
            DO itet = 1,kpts%ntet
               IF (ALL(kpts%ntetra(:,itet).ne.ikpt)) CYCLE

               eval = ev(iBand,kpts%ntetra(:,itet))

               IF(ANY(eval.GE.9999.9)) CYCLE

               DO i=1,4
                  icorn = kpts%ntetra(i,itet)

                  weight(i)=1.0
                  DO j=1,4
                     IF (i.NE.j) weight(i)=weight(i)*(eval(j)-eval(i))
                  ENDDO
                  weight(i)=6.0*kpts%voltet(itet)/(weight(i)*kpts%ntet)

                  DO idim=1,qdim
                     wpar(idim,iBand,icorn) =  wpar(idim,iBand,icorn) &
                                             + 0.25*weight(i)*qal(idim,iBand,ikpt)
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

      DO ikpt = 1,kpts%nkpt
         DO iBand = 1,neigd

            ener = ev(iBand,ikpt)
            DO idim = 1, qdim
               w  = 0.5*wpar(idim,iBand,ikpt)
               DO ie = 1,ned
                  term(ie) = energy(ie) - ener
                  IF(energy(ie).GT.ecmax(iBand)) term(ie) = ecmax(iBand) - ener
                  IF(term(ie).LT.0.0e0)         term(ie) = 0.0e0
                  g(ie,idim) = g(ie,idim) + w * term(ie)**2
               ENDDO
            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE tetra_dos
END MODULE m_tetrados
