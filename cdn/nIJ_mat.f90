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

    SUBROUTINE nIJ_mat(input,atoms,ne,usdus,jspin,we,eigVecCoeffs,cell,kpts,kptindx,nIJ_llp_mmp,enpara,v)

      USE m_types
      USE m_constants
      USE m_juDFT
      USE m_intgr, ONLY : intgr0
      USE m_radfun
      !USE m_check_mt_radii

      IMPLICIT NONE

      TYPE(t_usdus),       INTENT(IN)     :: usdus
      TYPE(t_input),       INTENT(IN)     :: input
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_eigVecCoeffs),INTENT(IN)     :: eigVecCoeffs
      TYPE(t_kpts),        INTENT(IN)     :: kpts
      TYPE(t_cell),        INTENT(IN)     :: cell
      INTEGER,             INTENT(IN)     :: ne,jspin,kptindx
      REAL,                INTENT(IN)     :: we(:)
      COMPLEX,             INTENT(INOUT)  :: nIJ_llp_mmp(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_enpara),      INTENT(IN)     :: enpara
      TYPE(t_potden),      INTENT(IN)     :: v
      TYPE(t_usdus)                       :: usdustemp 

      INTEGER i,i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2,counter
      COMPLEX c_0,A1, B1, A2, B2, power_factor, exponent
      REAL  norm1_W, norm2_W

        counter=0
        DO i_v = 1,atoms%n_v
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                counter=counter+1
            ENDDO
        ENDDO

        CALL usdustemp%init(atoms,input%jspins)
        CALL timestart("nIJ_mat")

        i_pair=1 
        DO i_v = 1,atoms%n_v
            natom1=atoms%lda_v(i_v)%atomIndex
            latom1=atoms%lda_v(i_v)%thisAtomL
            ll1atom1=latom1*(latom1+1)
            norm1_W = usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
                latom2=atoms%lda_v(i_v)%otherAtomL
                ll1atom2=latom2*(latom2+1)
                norm2_W = usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5
                power_factor=(cmplx(0, 1)**latom1) *(cmplx(0, -1)**latom2)
                !! power_factor is not included in the representation of matching coefficients in hsmt_ab.f90 routine. 
                !! Note that the $e^{ik.r_{atom I/J}}$ is included in c_ph(k,igSpin) only the shift exponent is needed.
                exponent=EXP(cmplx(0.0,-tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),kpts%bk(:,kptindx)))
                Do matom1=-latom1,latom1
                    lm1atom1=ll1atom1+matom1
                    Do matom2=-latom2,latom2
                        lm1atom2=ll1atom2+matom2
                        c_0=cmplx_0
                        Do i=1,ne 
                             A1 = eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)
                             B1 = eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)
                             A2 = eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin)
                             B2 = eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin)   
                             c_0 = c_0 + we(i) * (conjg(A2)*A1 + conjg(A2)*B1*norm1_W + conjg(B2)*A1*norm2_W + conjg(B2)*B1*norm1_W*norm2_W) *power_factor * exponent
                        ENDDO
                        nIJ_llp_mmp(matom1,matom2,i_pair) = nIJ_llp_mmp(matom1,matom2,i_pair) + c_0
                    ENDDO
                ENDDO
                i_pair=i_pair+1
            ENDDO
        ENDDO
        call timestop("nIJ_mat")
        
    END SUBROUTINE nIJ_mat
END MODULE m_nIJmat


