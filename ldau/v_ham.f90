MODULE m_vham

    !-------------------------------------------------------------------------------------------------------------|
    !     s                 s         ---   ---   IJ   IJ     s     J     I     s                                 |
    !<Phsi  |V_POT_MAT| Phis  > = -   >     >    V    n   <Psi  |phi > <phi |Psi  >                               |
    !     k                 k         ---   ---               k     Lp     L    k                                 |
    !                                 I,J,s  L,Lp                                                                 |
    !                                                                                                             |
    ! where L={l,m},Lp={lp,mp}                                                                                    |
    !                                                                                                             |
    !-------------------------------------------------------------------------------------------------------------|


    CONTAINS

    SUBROUTINE v_ham(input,usdus,atoms,kpts,cell,lapw,sym,noco,nococonv,fjgj,den,jspin,kptindx,hmat)

        USE m_types
        USE m_constants
        USE m_juDFT
        USE m_hsmt_ab
        USE m_hsmt_fjgj
        USE m_ylm
        USE m_radsrd


        IMPLICIT NONE

        TYPE(t_input),       INTENT(IN)     :: input
        TYPE(t_usdus),       INTENT(IN)     :: usdus
        TYPE(t_atoms),       INTENT(IN)     :: atoms
        TYPE(t_kpts),        INTENT(IN)     :: kpts
        TYPE(t_cell),        INTENT(IN)     :: cell
        TYPE(t_lapw),        INTENT(IN)     :: lapw
        TYPE(t_sym),         INTENT(IN)     :: sym
        TYPE(t_noco),        INTENT(IN)     :: noco
        TYPE(t_nococonv),    INTENT(IN)     :: nococonv
        TYPE(t_fjgj),        INTENT(INOUT)  :: fjgj
        TYPE(t_potden),      INTENT(IN)     :: den
        INTEGER,             INTENT(IN)     :: jspin,kptindx    
        CLASS(t_mat),        INTENT(INOUT)  :: hmat      !!G1, G2

        INTEGER i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2,iG1,iG2,abSizeG1,abSizeG2,indx_pair
        COMPLEX c_0,c_pair,c_hermitian
        COMPLEX, ALLOCATABLE :: abG1(:,:),abG2(:,:)
        COMPLEX, ALLOCATABLE :: c_0_pair(:,:,:)  !G1,G2,i_pair

        ALLOCATE(abG1(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
        ALLOCATE(abG2(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
        ALLOCATE(c_0_pair(atoms%lda_v(1)%numOtherAtoms * atoms%n_v ,MAXVAL(lapw%nv),MAXVAL(lapw%nv)))

        
        i_pair=1 
        DO i_v = 1,atoms%n_v  
            natom1=atoms%lda_v(i_v)%atomIndex
            latom1=atoms%lda_v(i_v)%thisAtomL
            ll1atom1=latom1*(latom1+1)
            CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom1),jspin)
            CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,jspin,atoms%itype(natom1),natom1,cell,lapw,fjgj,abG1,abSizeG1,.FALSE.)
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
                latom2=atoms%lda_v(i_v)%otherAtomL
                ll1atom2=latom2*(latom2+1)
                CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom2),jspin)
                CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,jspin,atoms%itype(natom2),natom2,cell,lapw,fjgj,abG2,abSizeG2,.FALSE.)
                DO iG1=1,lapw%nv(jspin)
                    Do iG2=1,lapw%nv(jspin)
                        c_0=cmplx_0
                        Do matom1=-latom1,latom1  
                            lm1atom1=ll1atom1+matom1
                            Do matom2=-latom2,latom2
                                lm1atom2=ll1atom2+matom2
                                !check pair density and which pair is sumed over take care of this 
                                c_0=c_0 - (atoms%lda_v(i_v)%V)*(conjg(den%nIJ_llp_mmp(matom1,matom2,i_pair,jspin)))*(conjg(abG1(lm1atom1+1,iG1))*abG2(lm1atom2+1,iG2) &
                                + conjg(abG1(lm1atom1+1+abSizeG1/2,iG1))*abG2(lm1atom2+1,iG2)*(usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5) &
                                + conjg(abG1(lm1atom1+1,iG1))*abG2(lm1atom2+1+abSizeG2/2,iG2)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5) &
                                + conjg(abG1(lm1atom1+1+abSizeG1/2,iG1))*abG2(lm1atom2+1+abSizeG2/2,iG2)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5)&
                                *(usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5))*&
                                EXP(cmplx(0.0,tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),(kpts%bk(:,kptindx)+lapw%gvec(:, iG2,jspin)))) &
                                *(cmplx(0, -1)**latom1) *(cmplx(0, 1)**latom2) 
                            ENDDO
                        ENDDO
                        c_0_pair(i_pair,iG1,iG2)=c_0 
                    ENDDO
                ENDDO
                i_pair=i_pair+1
            ENDDO
        ENDDO

        DO iG1=1,lapw%nv(jspin)
            Do iG2=1,lapw%nv(jspin)
                c_pair=cmplx_0
                c_hermitian=cmplx_0
                DO indx_pair=1, i_pair-1
                    c_pair=c_pair+c_0_pair(indx_pair,iG1,iG2)
                    c_hermitian=c_hermitian+conjg(c_0_pair(indx_pair,iG2,iG1))
                ENDDO
                WRITE(6000,*) 'iG1,iG2,mat',iG1,iG2,c_pair
                IF (iG1==iG2) THEN
                    WRITE(3000,*) 'iG1,iG2,mat',iG1,iG2,c_pair
                ENDIF
                WRITE(3500,*) 'iG2,iG1,difference', iG2,iG1,c_hermitian - c_pair
                IF(hmat%l_real) THEN
                    hmat%data_r(iG1,iG2)=hmat%data_r(iG1,iG2)+REAL(c_pair)
                ELSE
                    hmat%data_c(iG1,iG2)=hmat%data_c(iG1,iG2)+c_pair
                END IF
            ENDDO
        ENDDO
    END SUBROUTINE v_ham
END MODULE m_vham