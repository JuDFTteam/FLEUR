
MODULE m_hyb_abcrot
CONTAINS
   SUBROUTINE hyb_abcrot(hybrid, atoms, neig, sym, cell, oneD,&
  &                 acof, bcof, ccof)
!     ***************************************************************
!     * This routine transforms a/b/cof which are given wrt rotated *
!     * MT functions (according to invsat/ngopr) into a/b/cof wrt   *
!     * unrotated MT functions. Needed for GW calculations.         *
!     *                                                             *
!     * Christoph Friedrich Mar/2005                                *
!     ***************************************************************
      USE m_dwigner
      USE m_types
      IMPLICIT NONE
      TYPE(t_hybrid), INTENT(IN) :: hybrid
      TYPE(t_oneD), INTENT(IN)   :: oneD
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_atoms), INTENT(IN)  :: atoms
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: neig
!     ..
!     .. Array Arguments ..

      COMPLEX, INTENT(INOUT) :: acof(:, 0:, :) !(dimension%neigd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT(INOUT) :: bcof(:, 0:, :) !(dimension%neigd,0:dimension%lmd,atoms%natd)
      COMPLEX, INTENT(INOUT) :: ccof(-atoms%llod:, :, :, :)!(-llod:llod,dimension%neigd,atoms%nlod,atoms%natd)
!     ..
!     .. Local Scalars ..
      INTEGER itype, ineq, iatom, iop, ilo, i, l, lm, lmp, ifac
!     ..
!     .. Local Arrays ..
!***** COMPLEX, ALLOCATABLE :: d_wgn(:,:,:,:) !put into module m_savewigner
!

      IF (.NOT. ALLOCATED(hybrid%d_wgn2)) THEN    !calculate sym%d_wgn only once
#ifndef CPP_MPI
         PRINT *, "calculate wigner-matrix"
#endif
         STOP "WIGNER MATRIX should be available in hybrid part"
         !IF (.NOT.oneD%odi%d1) THEN
         !  ALLOCATE (sym%d_wgn(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,sym%nop))
         !  CALL d_wigner(sym%nop,sym%mrot,cell%bmat,atoms%lmaxd,sym%d_wgn)
         !ELSE
         !  ALLOCATE (sym%d_wgn(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,atoms%lmaxd,oneD%ods%nop))
         !  CALL d_wigner(oneD%ods%nop,oneD%ods%mrot,cell%bmat,atoms%lmaxd,sym%d_wgn)
         !ENDIF
      ENDIF

      iatom = 0
      DO itype = 1, atoms%ntype
         DO ineq = 1, atoms%neq(itype)
            iatom = iatom + 1
            IF (.NOT. oneD%odi%d1) THEN
               iop = atoms%ngopr(iatom)
            ELSE
               iop = oneD%ods%ngopr(iatom)
            ENDIF
!                                    l                        l    l
! inversion of spherical harmonics: Y (pi-theta,pi+phi) = (-1)  * Y (theta,phi)
!                                    m                             m
            ifac = 1
            IF (atoms%invsat(iatom) == 2) THEN
               IF (.NOT. oneD%odi%d1) THEN
                  iop = atoms%ngopr(sym%invsatnr(iatom))
               ELSE
                  iop = oneD%ods%ngopr(sym%invsatnr(iatom))
               ENDIF
               ifac = -1
            ENDIF
            DO l = 1, atoms%lmax(itype)
!  replaced d_wgn by conjg(d_wgn),FF October 2006
               DO i = 1, neig
                  acof(i, l**2:l*(l + 2), iatom) = ifac**l*matmul(conjg(hybrid%d_wgn2(-l:l, -l:l, l, iop)), acof(i, l**2:l*(l + 2), iatom))
                  bcof(i, l**2:l*(l + 2), iatom) = ifac**l*matmul(conjg(hybrid%d_wgn2(-l:l, -l:l, l, iop)), bcof(i, l**2:l*(l + 2), iatom))
               ENDDO
            ENDDO
            DO ilo = 1, atoms%nlo(itype)
               l = atoms%llo(ilo, itype)
               IF (l > 0) THEN
                  DO i = 1, neig
                     ccof(-l:l, i, ilo, iatom) = ifac**l*matmul(conjg(hybrid%d_wgn2(-l:l, -l:l, l, iop)), ccof(-l:l, i, ilo, iatom))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE hyb_abcrot
END MODULE m_hyb_abcrot
