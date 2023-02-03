MODULE m_nIJmat
     
    CONTAINS

    SUBROUTINE nIJ_mat(atoms,ne,usdus,jspin,we,eigVecCoeffs,cell,kpts,kptindx,nIJ_llp_mmp)

      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_usdus),       INTENT(IN)     :: usdus
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_eigVecCoeffs),INTENT(IN)     :: eigVecCoeffs
      TYPE(t_cell),        INTENT(IN)     :: cell
      TYPE(t_kpts),        INTENT(IN)     :: kpts
      INTEGER,             INTENT(IN)     :: ne,jspin,kptindx
      REAL,                INTENT(IN)     :: we(:)
      COMPLEX,             INTENT(INOUT)  :: nIJ_llp_mmp(:,:,:)!!(-latom1:latom1(to save m1),-latom2:latom2,ii_v)

      INTEGER i,i_v,ii_v,atom1,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2
      REAL shft1,shft2,shft3
      COMPLEX c_0

      ii_v=1 !counts number of pairs
        DO i_v = 1,atoms%n_v  !loop over pairs which are corrected by U+V 
            natom1=atoms%lda_v(i_v)%atomIndex
            latom1=atoms%lda_v(i_v)%thisAtomL
            ll1atom1=latom1*(latom1+1)
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms,1
                natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
                latom2=atoms%lda_v(i_v)%otherAtomL
                ll1atom2=latom2*(latom2+1)
                shft1=atoms%lda_v(i_v)%atomShifts(1,atom2) * kpts%bk(1,kptindx) * (cell%amat(1,1)**2 + cell%amat(2,1)**2 + cell%amat(3,1)**2)**0.5
                shft2=atoms%lda_v(i_v)%atomShifts(2,atom2) * kpts%bk(2,kptindx) * (cell%amat(1,2)**2 + cell%amat(2,2)**2 + cell%amat(3,2)**2)**0.5
                shft3=atoms%lda_v(i_v)%atomShifts(3,atom2) * kpts%bk(3,kptindx) * (cell%amat(1,3)**2 + cell%amat(2,3)**2 + cell%amat(3,3)**2)**0.5
                Do matom1=-latom1,latom1,1
                    lm1atom1=ll1atom1+matom1
                    Do matom2=-latom2,latom2,1
                        lm1atom2=ll1atom2+matom2
                        c_0=cmplx_0
                        Do i=1,ne    
                             c_0 = c_0 + we(i) * ( conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom1,natom1,jspin)**0.5) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)*(usdus%ddn(latom2,natom2,jspin)**0.5) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom2,natom2,jspin)**0.5) &
                             *(usdus%ddn(latom1,natom1,jspin))**0.5) * EXP(cmplx(0.0,-1.0)*(shft1+shft2+shft3))
                        ENDDO
                        nIJ_llp_mmp(matom1,matom2,ii_v)= c_0
                    ENDDO
                ENDDO
                ii_v=ii_v+1
            ENDDO
        ENDDO
    END SUBROUTINE nIJ_mat
END MODULE m_nIJmat