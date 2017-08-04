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
   INTEGER                 :: nsym,ID_mat(3,3)
   
   nsym=sym%nop
   if (.not.sym%invs) nsym=2*sym%nop
   
   ALLOCATE (kpts%bkf(3,nsym*kpts%nkpt))
   ALLOCATE (kpts%bkp(nsym*kpts%nkpt))
   ALLOCATE (kpts%bksym(nsym*kpts%nkpt))
      
   ! Generate symmetry operations in reciprocal space
   DO iop=1,nsym
      IF( iop .le. sym%nop ) THEN
         rrot(:,:,iop) = transpose( sym%mrot(:,:,sym%invtab(iop)) )
      ELSE
         rrot(:,:,iop) = -rrot(:,:,iop-sym%nop)
      END IF
   END DO
  
   !Add existing vectors to list of full vectors
   id_mat=0
   ID_mat(1,1)=1;ID_mat(2,2)=1;ID_mat(3,3)=1
   IF (ANY(sym%mrot(:,:,1).NE.ID_mat)) CALL judft_error("Identity must be first symmetry operation",calledby="gen_bz")
   
   ic=0
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
