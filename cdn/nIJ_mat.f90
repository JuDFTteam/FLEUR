MODULE m_nIJmat
     
    CONTAINS

    SUBROUTINE nIJ_mat(atoms,ne,usdus,jspin,we,eigVecCoeffs,cell,kpts,kptindx,nIJ_llp_mmp)

      USE m_types
      USE m_constants
      USE m_juDFT

      IMPLICIT NONE

      TYPE(t_usdus),       INTENT(IN)     :: usdus
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_eigVecCoeffs),INTENT(IN)     :: eigVecCoeffs
      TYPE(t_kpts),        INTENT(IN)     :: kpts
      TYPE(t_cell),        INTENT(IN)     :: cell
      INTEGER,             INTENT(IN)     :: ne,jspin,kptindx
      REAL,                INTENT(IN)     :: we(:)
      COMPLEX,             INTENT(INOUT)  :: nIJ_llp_mmp(-lmaxU_const:,-lmaxU_const:,:)
     

      INTEGER i,i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2
      COMPLEX c_0
      COMPLEX nII_ll_mmp

        call timestart("nIJ_mat")
        i_pair=1 !counts number of pairs
        DO i_v = 1,atoms%n_v!1  
            natom1=atoms%lda_v(i_v)%atomIndex!1
            latom1=atoms%lda_v(i_v)%thisAtomL!1
            ll1atom1=latom1*(latom1+1)
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)!1
                latom2=atoms%lda_v(i_v)%otherAtomL!1
                ll1atom2=latom2*(latom2+1)
                Do matom1=-latom1,latom1
                    lm1atom1=ll1atom1+matom1
                    Do matom2=-latom2,latom2
                        lm1atom2=ll1atom2+matom2
                        c_0=cmplx_0
                        Do i=1,ne    
                             c_0 = c_0 + we(i) * (conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)&
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5) &
                             *(usdus%ddn(latom1,atoms%itype(natom1),jspin))**0.5)* EXP(cmplx(0.0,-tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),kpts%bk(:,kptindx)))
                             WRITE(2000,*) 'aa',conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)
                             WRITE(2000,*) 'ab',conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5)
                             WRITE(2000,*) 'ba',conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5)
                             WRITE(2000,*) 'bb',conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5)*(usdus%ddn(latom1,atoms%itype(natom1),jspin))**0.5
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
