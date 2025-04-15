MODULE m_nIJmat

   !------------------------------------------------------------------------------ !
   !                                                                               !
   ! MODULE: m_nIJmat                                                              !
   !                                                                               !
   ! Author : Wejdan Beida                                                         !
   !                                                                               !
   ! Description: This module calculates the occupation matrix between two atoms.  !
   !                                                                               !
   !------------------------------------------------------------------------------ !

CONTAINS

   SUBROUTINE nIJ_mat(firstspin, input, atoms, ne, usdus, we, abc, cell, kpts, kptindx, nIJ_llp_mmp, enpara, v)

      USE m_types
      USE m_constants
      USE m_juDFT
      USE m_intgr, ONLY: intgr0
      USE m_radfun
      USE m_types_abc
      !USE m_check_mt_radii

      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: firstspin
      TYPE(t_usdus), INTENT(IN)    :: usdus
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_abc), INTENT(IN)      :: abc(:, firstspin:)
      TYPE(t_kpts), INTENT(IN)     :: kpts
      TYPE(t_cell), INTENT(IN)     :: cell
      INTEGER, INTENT(IN)          :: ne,  kptindx
      REAL, INTENT(IN)             :: we(:)
      COMPLEX, INTENT(INOUT)       :: nIJ_llp_mmp(-lmaxU_const:, -lmaxU_const:, :, :)
      TYPE(t_enpara), INTENT(IN)   :: enpara
      TYPE(t_potden), INTENT(IN)   :: v
      TYPE(t_usdus)                :: usdustemp

      INTEGER i,i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2,counter,itype1,itype2,jspin
      COMPLEX c_0, A1, B1, A2, B2, power_factor, exponent
      REAL norm1_W, norm2_W, r_5digit, i_5digit
      counter = 0
      DO i_v = 1, atoms%n_v
         Do atom2 = 1, atoms%lda_v(i_v)%numOtherAtoms
            counter = counter + 1
         END DO
      END DO

      CALL usdustemp%init(atoms, input%jspins)
      CALL timestart("nIJ_mat")

      DO jspin = lbound(abc, 2), ubound(abc, 2)
         i_pair = 1
         DO i_v = 1, atoms%n_v
            natom1 = atoms%lda_v(i_v)%atomIndex
            itype1 = atoms%itype(natom1)
            natom1 = natom1 - atoms%firstatom(itype1) + 1
            latom1 = atoms%lda_v(i_v)%thisAtomL
            ll1atom1 = latom1*(latom1 + 1)
            norm1_W = usdus%ddn(latom1, atoms%itype(natom1), jspin)**0.5
            Do atom2 = 1, atoms%lda_v(i_v)%numOtherAtoms
               natom2 = atoms%lda_v(i_v)%otherAtomIndices(atom2)
               itype2 = atoms%itype(natom2)
               natom2 = natom2 - atoms%firstatom(itype2) + 1
               latom2 = atoms%lda_v(i_v)%otherAtomL
               ll1atom2 = latom2*(latom2 + 1)
               norm2_W = usdus%ddn(latom2, atoms%itype(natom2), jspin)**0.5
               power_factor = (cmplx(0, 1)**latom1)*(cmplx(0, -1)**latom2)
                !! power_factor is not included in the representation of matching coefficients in hsmt_ab.f90 routine.
                !! Note that the $e^{ik.r_{atom I/J}}$ is included in c_ph(k,igSpin) only the shift exponent is needed.
               exponent = EXP(cmplx(0.0, -tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:, atom2), kpts%bk(:, kptindx)))
               Do matom1 = -latom1, latom1
                  lm1atom1 = ll1atom1 + matom1
                  Do matom2 = -latom2, latom2
                     lm1atom2 = ll1atom2 + matom2
                     c_0 = cmplx_0
                     Do i = 1, ne
                        A1 = abc(itype1, jspin)%cof(i, lm1atom1, 1, natom1)
                        B1 = abc(itype1, jspin)%cof(i, lm1atom1, 2, natom1)
                        A2 = abc(itype2, jspin)%cof(i, lm1atom2, 1, natom2)
                        B2 = abc(itype2, jspin)%cof(i, lm1atom2, 2, natom2)
                        c_0 = c_0 + we(i)*(conjg(A2)*A1 + conjg(A2)*B1*norm1_W + &
                                           conjg(B2)*A1*norm2_W + conjg(B2)*B1*norm1_W*norm2_W)*power_factor*exponent
                     END DO
                     nIJ_llp_mmp(matom1, matom2, i_pair, jspin) = nIJ_llp_mmp(matom1, matom2, i_pair,jspin) + c_0
                     !r_5digit= anint(REAL(nIJ_llp_mmp(matom1,matom2,i_pair))*10.0**5)/10.0**5
                     !i_5digit=anint(AIMAG(nIJ_llp_mmp(matom1,matom2,i_pair))*10.0**5)/10.0**5
                     !nIJ_llp_mmp(matom1,matom2,i_pair) = cmplx(r_5digit,i_5digit)
                     !WRITE(211,*) "nIJ", nIJ_llp_mmp(matom1,matom2,i_pair)
                  END DO
               END DO
               i_pair = i_pair + 1
            END DO
         END DO
      end do
      call timestop("nIJ_mat")

   END SUBROUTINE nIJ_mat
END MODULE m_nIJmat

