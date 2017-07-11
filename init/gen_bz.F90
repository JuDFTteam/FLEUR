!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gen_bz generates the (whole) Brillouin zone from the          !
! (irreducible) k-points given in the kpts file.                !
!                                                               !
!                                     M.Betzinger (09/07)       !
!                                                               !
!                        Refactored in 2017 by G.M.             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
MODULE m_gen_bz

CONTAINS

SUBROUTINE gen_bz( kpts,sym)

   !     bk     ::    irreducible k-points
   !     nkpt   ::    number of irr. k-points
   !     bkf    ::    all k-points
   !     nkptf  ::    number of all k-points
   !     bkp    ::    k-point parent
   !     bksym  ::    symmetry operation, that connects the parent
   !                  k-point with the current one

   USE m_juDFT
   USE m_util, ONLY: modulo1
   USE m_types

   IMPLICIT NONE

   TYPE(t_kpts),INTENT(INOUT) :: kpts
   TYPE(t_sym),INTENT(IN)     :: sym

!  - local scalars -
   INTEGER                 ::  ic,iop,ikpt,ikpt1
   LOGICAL                 ::  l_found
      
!  - local arrays - 
   INTEGER,ALLOCATABLE     ::  iarr(:)
   REAL                    ::  rrot(3,3,2*sym%nop),rotkpt(3)
   REAL,ALLOCATABLE        ::  rarr1(:,:)
   INTEGER:: nsym
   
   nsym=sym%nop
   if (.not.sym%invs) nsym=2*sym%nop
   
   ALLOCATE (kpts%bkf(3,nsym*kpts%nkpt))
   ALLOCATE (kpts%bkp(nsym*kpts%nkpt))
   ALLOCATE (kpts%bksym(nsym*kpts%nkpt))
      
   ! Generate symmetry operations in reciprocal space
   DO iop=1,nsym
      IF( iop .le. sym%nop ) THEN
         rrot(:,:,iop) = transpose( sym%mrot(:,:,iop) )
      ELSE
         rrot(:,:,iop) = -rrot(:,:,iop-sym%nop)
      END IF
   END DO

   ! Set target number for k points in full BZ
   kpts%nkptf = kpts%nkpt3(1)*kpts%nkpt3(2)*kpts%nkpt3(3)
   IF(kpts%l_gamma) THEN
      IF (ANY(MODULO(kpts%nkpt3(:),2).EQ.0)) THEN
         kpts%nkptf = kpts%nkptf + 1
      END IF
   END IF

   ! Apply symmetrie operations to all k-points of IBZ, test whether
   ! generated k-point already is in the full BZ set of k-points, and
   ! add it if it is not yet in this set.

   kpts%bkf = 0
  
   !Add existing vectors to list of full vectors
   print *,"WARNING from gen_bz"
   print *,"Assuming Identity to be fist symmetry op!"
   DO ic=1,kpts%nkpt
      kpts%bkf(:,ic) = kpts%bk(:,ic)
      kpts%bkp(ic)  = ic
      kpts%bksym(ic) = 1
   ENDDO
   ic=ic-1
   
   DO iop=1,nsym
      DO ikpt=1,kpts%nkpt
         l_found = .FALSE.
         rotkpt = MATMUL(rrot(:,:,iop), kpts%bk(:,ikpt))
         !transform back into IBZ
         rotkpt = modulo1(rotkpt,kpts%nkpt3)
         DO ikpt1=1,ic
            IF (MAXVAL(ABS(kpts%bkf(:,ikpt1) - rotkpt)).LE.1e-08) THEN
               l_found = .TRUE.
               EXIT
            END IF
         END DO
          
         IF(.NOT.l_found) THEN
            ic = ic + 1
            kpts%bkf(:,ic) = rotkpt
            kpts%bkp(ic) = ikpt
            kpts%bksym(ic) = iop
         END IF
      END DO
   END DO

   kpts%nkptf = ic

   ! Reallocate bkf, bkp, bksym
   ALLOCATE (iarr(kpts%nkptf))
   iarr = kpts%bkp(:kpts%nkptf)
   DEALLOCATE(kpts%bkp)
   ALLOCATE (kpts%bkp(kpts%nkptf))
   kpts%bkp = iarr
   iarr= kpts%bksym(:kpts%nkptf)
   DEALLOCATE (kpts%bksym )
   ALLOCATE (kpts%bksym(kpts%nkptf))
   kpts%bksym = iarr
   DEALLOCATE(iarr)
   ALLOCATE (rarr1(3,kpts%nkptf))
   rarr1 = kpts%bkf(:,:kpts%nkptf)
   DEALLOCATE (kpts%bkf )
   ALLOCATE (kpts%bkf(3,kpts%nkptf))
   kpts%bkf = rarr1
   DEALLOCATE(rarr1)
      
END SUBROUTINE gen_bz

END MODULE m_gen_bz
