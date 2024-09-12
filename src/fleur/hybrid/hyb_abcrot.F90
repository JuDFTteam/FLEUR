MODULE m_hyb_abcrot
CONTAINS
   SUBROUTINE hyb_abcrot(hybinp, atoms, neig, sym,&
                    acof, bcof, ccof)
!     ***************************************************************
!     * This routine transforms a/b/cof which are given wrt rotated *
!     * MT functions (according to invsat/ngopr) into a/b/cof wrt   *
!     * unrotated MT functions. Needed for GW calculations.         *
!     *                                                             *
!     * Christoph Friedrich Mar/2005                                *
!     ***************************************************************
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_atoms), INTENT(IN)  :: atoms
      INTEGER, INTENT(IN) :: neig

      COMPLEX, INTENT(INOUT) :: acof(:, 0:, :) !(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd)
      COMPLEX, INTENT(INOUT) :: bcof(:, 0:, :) !(input%neig,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd)
      COMPLEX, INTENT(INOUT) :: ccof(-atoms%llod:, :, :, :)!(-llod:llod,input%neig,atoms%nlod,atoms%natd)
      INTEGER itype, ineq, iatom, iop, ilo, i, l, ifac
      call timestart("hyb_abcrot")
      IF (.NOT. ALLOCATED(hybinp%d_wgn2)) THEN    !calculate sym%d_wgn only once
         PRINT *, "calculate wigner-matrix"
         call judft_error('WIGNER MATRIX should be available in hybinp part')
      ENDIF

      !$OMP PARALLEL DO default(none) private(iatom, itype, iop, ifac, l, i, ilo) &
      !$OMP shared(atoms, sym, acof, bcof, ccof, neig, hybinp)
      do iatom = 1, atoms%nat 
         itype = atoms%itype(iatom)
         iop = sym%ngopr(iatom)
         !                                    l                        l    l
         ! inversion of spherical harmonics: Y (pi-theta,pi+phi) = (-1)  * Y (theta,phi)
         !                                    m                             m
         ifac = 1
         IF (sym%invsat(iatom) == 2) THEN
            iop = sym%ngopr(sym%invsatnr(iatom))

            ifac = -1
         ENDIF
         DO l = 1, atoms%lmax(itype)
            DO i = 1, neig
               acof(i, l**2:l*(l + 2), iatom) = ifac**l*matmul(conjg(hybinp%d_wgn2(-l:l, -l:l, l, iop)),acof(i, l**2:l*(l + 2), iatom))
               bcof(i, l**2:l*(l + 2), iatom) = ifac**l*matmul(conjg(hybinp%d_wgn2(-l:l, -l:l, l, iop)),bcof(i, l**2:l*(l + 2), iatom))
            ENDDO
         ENDDO
         DO ilo = 1, atoms%nlo(itype)
            l = atoms%llo(ilo, itype)
            IF (l > 0) THEN
               DO i = 1, neig
                  ccof(-l:l, i, ilo, iatom) = ifac**l*matmul(conjg(hybinp%d_wgn2(-l:l, -l:l, l, iop)), ccof(-l:l, i, ilo, iatom))
               ENDDO
            ENDIF
            ENDDO
      ENDDO
      !$OMP end parallel do
      call timestop("hyb_abcrot")
   END SUBROUTINE hyb_abcrot
END MODULE m_hyb_abcrot
