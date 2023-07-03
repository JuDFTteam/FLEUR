MODULE m_nIJmat
     
    CONTAINS

    SUBROUTINE nIJ_mat(input,atoms,ne,usdus,jspin,we,eigVecCoeffs,cell,kpts,kptindx,nIJ_llp_mmp,enpara,v)

      USE m_types
      USE m_constants
      USE m_juDFT
      USE m_intgr, ONLY : intgr0
      USE m_radfun

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

      REAL, ALLOCATABLE :: f1(:,:), g1(:,:),u1rsquare_l(:),udotrsquare_l(:),u2rsquare_l(:),udotrsquare2_l(:),f2(:,:), g2(:,:),u1rsquare_s(:),udotrsquare_s(:),u2rsquare_s(:),udotrsquare2_s(:)
      COMPLEX, ALLOCATABLE :: henning_mat(:,:,:)
      REAL intu1_l,intudot1_l,wronk1,wronk2,intu2_l,intudot2_l,intu1_s,intudot1_s,intu2_s,intudot2_s
      INTEGER j,nodeu1,noded1,nodeu2,noded2


      INTEGER i,i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2
      COMPLEX c_0,h_0

      ALLOCATE(f1(atoms%jmtd,2),g1(atoms%jmtd,2),f2(atoms%jmtd,2),g2(atoms%jmtd,2))
      ALLOCATE(u1rsquare_l(atoms%jmtd),udotrsquare_l(atoms%jmtd),u2rsquare_l(atoms%jmtd),udotrsquare2_l(atoms%jmtd))
      ALLOCATE(u1rsquare_s(atoms%jmtd),udotrsquare_s(atoms%jmtd),u2rsquare_s(atoms%jmtd),udotrsquare2_s(atoms%jmtd))
      ALLOCATE(henning_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_v))

        CALL usdustemp%init(atoms,input%jspins)
        CALL timestart("nIJ_mat")
        i_pair=1 !counts number of pairs
        DO i_v = 1,atoms%n_v!1  
            natom1=atoms%lda_v(i_v)%atomIndex!1
            latom1=atoms%lda_v(i_v)%thisAtomL!1
            ll1atom1=latom1*(latom1+1)
            !!!
            CALL radfun(latom1,atoms%itype(natom1),jspin,enpara%el0(latom1,atoms%itype(natom1),jspin),v%mt(:,0,atoms%itype(natom1),jspin),atoms,f1(:,:),g1(:,:),usdustemp,nodeu1,noded1,wronk1)
            DO j=1,atoms%jri(atoms%itype(natom1))
                u1rsquare_l(j)   = f1(j,1)*atoms%rmsh(j,atoms%itype(natom1))
                udotrsquare_l(j) = g1(j,1)*atoms%rmsh(j,atoms%itype(natom1))
                u1rsquare_s(j)   = f1(j,2)*atoms%rmsh(j,atoms%itype(natom1))
                udotrsquare_s(j) = g1(j,2)*atoms%rmsh(j,atoms%itype(natom1))
            ENDDO
            CALL intgr0(u1rsquare_l,atoms%rmsh(1,atoms%itype(natom1)),atoms%dx(atoms%itype(natom1)),atoms%jri(atoms%itype(natom1)),intu1_l)
            CALL intgr0(udotrsquare_l,atoms%rmsh(1,atoms%itype(natom1)),atoms%dx(atoms%itype(natom1)),atoms%jri(atoms%itype(natom1)),intudot1_l)
            CALL intgr0(u1rsquare_s,atoms%rmsh(1,atoms%itype(natom1)),atoms%dx(atoms%itype(natom1)),atoms%jri(atoms%itype(natom1)),intu1_s)
            CALL intgr0(udotrsquare_s,atoms%rmsh(1,atoms%itype(natom1)),atoms%dx(atoms%itype(natom1)),atoms%jri(atoms%itype(natom1)),intudot1_s)
            WRITE(66,*) 'int u1* r^2 dr , int u1dot* r^2 dr ', intu1_l,intu1_s, intudot1_l, intudot1_s 
            !!!
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)!1
                latom2=atoms%lda_v(i_v)%otherAtomL!1
                ll1atom2=latom2*(latom2+1)
                !!!
                CALL radfun(latom2,atoms%itype(natom2),jspin,enpara%el0(latom2,atoms%itype(natom2),jspin),v%mt(:,0,atoms%itype(natom2),jspin),atoms,f2(:,:),g2(:,:),usdustemp,nodeu2,noded2,wronk2)
                DO j = 1,atoms%jri(atoms%itype(natom2))
                    u2rsquare_l(j)    = f2(j,1)*atoms%rmsh(j,atoms%itype(natom2))
                    udotrsquare2_l(j) = g2(j,1)*atoms%rmsh(j,atoms%itype(natom2))
                    u2rsquare_s(j)    = f2(j,2)*atoms%rmsh(j,atoms%itype(natom2))
                    udotrsquare2_s(j) = g2(j,2)*atoms%rmsh(j,atoms%itype(natom2))
                END DO
                CALL intgr0(u2rsquare_l,atoms%rmsh(1,atoms%itype(natom2)),atoms%dx(atoms%itype(natom2)),atoms%jri(atoms%itype(natom2)),intu2_l)
                CALL intgr0(udotrsquare2_l,atoms%rmsh(1,atoms%itype(natom2)),atoms%dx(atoms%itype(natom2)),atoms%jri(atoms%itype(natom2)),intudot2_l)
                CALL intgr0(u2rsquare_s,atoms%rmsh(1,atoms%itype(natom2)),atoms%dx(atoms%itype(natom2)),atoms%jri(atoms%itype(natom2)),intu2_s)
                CALL intgr0(udotrsquare2_s,atoms%rmsh(1,atoms%itype(natom2)),atoms%dx(atoms%itype(natom2)),atoms%jri(atoms%itype(natom2)),intudot2_s)
                WRITE(66,*) 'int u2* r^2 dr , int u2dot* r^2 dr ', intu2_l, intu2_s, intudot2_l, intudot2_s
                !!!
                Do matom1=-latom1,latom1
                    lm1atom1=ll1atom1+matom1
                    Do matom2=-latom2,latom2
                        lm1atom2=ll1atom2+matom2
                        c_0=cmplx_0
                        h_0=cmplx_0
                        Do i=1,ne    
                             c_0 = c_0 + we(i) * (conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)&
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5) &
                             *(usdus%ddn(latom1,atoms%itype(natom1),jspin))**0.5)* EXP(cmplx(0.0,-tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),kpts%bk(:,kptindx)))&
                             *(cmplx(0, 1)**latom1) *(cmplx(0, -1)**latom2) 
                            !!!!!!
                            h_0= h_0 + we(i) * (conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)*((intu1_l*intu2_l)+(intu1_s*intu2_s))&
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,0,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*((intudot1_l*intu2_l)+(intudot1_s*intu2_s)) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,0,natom1,jspin)*((intudot2_l*intu1_l)+(intudot2_s*intu1_s)) &
                             +conjg(eigVecCoeffs%abcof(i,lm1atom2,1,natom2,jspin))*eigVecCoeffs%abcof(i,lm1atom1,1,natom1,jspin)*((intudot1_l*intudot2_l)+(intudot1_s*intudot2_s))) &
                             * EXP(cmplx(0.0,-tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),kpts%bk(:,kptindx)))&
                             *(cmplx(0, 1)**latom1) *(cmplx(0, -1)**latom2)
                            !!!!!!
                        ENDDO
                        nIJ_llp_mmp(matom1,matom2,i_pair) = nIJ_llp_mmp(matom1,matom2,i_pair) + c_0
                        henning_mat(matom1,matom2,i_pair) = henning_mat(matom1,matom2,i_pair) + h_0
                        WRITE(33,*) 'Wejdan_contrib, Henning_contrib', c_0, h_0
                        WRITE(444,*) 'pair,m1,m2,Heninng Mat', i_pair,matom1, matom2,henning_mat(matom1,matom2,i_pair)
                        WRITE(555,*) 'pair,m1,m2,Wejdan Mat', i_pair,matom1, matom2,nIJ_llp_mmp(matom1,matom2,i_pair)
                    ENDDO
                ENDDO
                WRITE(222,*) 'pair index, term1 Wejdan, term1 Henning', i_pair, 1, (intu2_l+intu2_s)*(intu1_l+intu1_s)
                WRITE(222,*) 'pair index, term2 Wejdan, term2 Henning', i_pair, usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5, (intu2_l+intu2_s)*(intudot1_l+intudot1_s)
                WRITE(222,*) 'pair index, term3 Wejdan, term3 Henning', i_pair, usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5 , (intu1_l+intu1_s)*(intudot2_l+intudot2_s)
                WRITE(222,*) 'pair index, term4 Wejdan, term4 Henning', i_pair, (usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5)*((usdus%ddn(latom1,atoms%itype(natom1),jspin))**0.5) ,(intudot1_l+intudot1_s)*(intudot2_l+intudot2_s)
                i_pair=i_pair+1
            ENDDO
        ENDDO
        call timestop("nIJ_mat")

    END SUBROUTINE nIJ_mat
END MODULE m_nIJmat

!!CALL radfun(latom1,atoms%itype(natom1),jspin,enpara%el0(latom1,atoms%itype(natom1),jspin),v%mt(:,0,atoms%itype(natom1),jspin),atoms,f1(:,2),g1(:,2),usdustemp,nodeu1_s,noded1_s,wronk1_s)
!!CALL radfun(latom2,atoms%itype(natom2),jspin,enpara%el0(latom2,atoms%itype(natom2),jspin),v%mt(:,0,atoms%itype(natom2),jspin),atoms,f2(:,1),g2(:,1),usdustemp,nodeu2_s,noded2_s,wronk2_s)