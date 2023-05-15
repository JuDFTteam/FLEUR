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

    SUBROUTINE v_ham(usdus,atoms,kpts,cell,lapw,sym,noco,nococonv,fjgj,den,jspin,kptindx,hmat)

        USE m_types
        USE m_constants
        USE m_juDFT
        USE m_hsmt_ab
        USE m_hsmt_fjgj
        USE m_ylm
        !USE m_matmul_dgemm


        IMPLICIT NONE

        TYPE(t_usdus),       INTENT(IN)     :: usdus
        TYPE(t_atoms),       INTENT(IN)     :: atoms
        TYPE(t_kpts),        INTENT(IN)     :: kpts
        TYPE(t_cell),        INTENT(IN)     :: cell
        TYPE(t_lapw),        INTENT(IN)     :: lapw
        TYPE(t_sym),         INTENT(IN)     :: sym
        TYPE(t_noco),        INTENT(IN)     :: noco
        TYPE(t_nococonv),    INTENT(IN)     :: nococonv
        TYPE(t_fjgj),        INTENT(IN)     :: fjgj
        TYPE(t_potden),      INTENT(IN)     :: den
        INTEGER,             INTENT(IN)     :: jspin,kptindx    
        CLASS(t_mat),        INTENT(INOUT)  :: hmat      !!G1, G2,i_pair,jspin

        INTEGER i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2,iG1,iG2,abSizeG1,abSizeG2
        COMPLEX c_0
        COMPLEX, ALLOCATABLE :: abG1(:,:),abG2(:,:)
        ALLOCATE(abG1(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
        ALLOCATE(abG2(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))

        i_pair=1 
        DO i_v = 1,atoms%n_v  
            natom1=atoms%lda_v(i_v)%atomIndex
            latom1=atoms%lda_v(i_v)%thisAtomL
            ll1atom1=latom1*(latom1+1)
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
                latom2=atoms%lda_v(i_v)%otherAtomL
                ll1atom2=latom2*(latom2+1)
                CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,jspin,natom1,natom1,cell,lapw,fjgj,abG1,abSizeG1,.FALSE.)
                CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,jspin,natom2,natom2,cell,lapw,fjgj,abG2,abSizeG2,.FALSE.)
                DO iG1=1,lapw%nv(jspin)
                    Do iG2=1,lapw%nv(jspin)
                        c_0=cmplx_0
                        Do matom1=-latom1,latom1  
                            lm1atom1=ll1atom1+matom1
                            Do matom2=-latom2,latom2
                                lm1atom2=ll1atom2+matom2
                                !!!!!!!!c_0=c_0 - atoms%lda_v(i_v)%V*(conjg(abG1(ll1atom1+1,iG1))*abG2(ll1atom2+1,iG2))! & !!!!!!!
                                c_0=c_0 - (atoms%lda_v(i_v)%V)*(den%nIJ_llp_mmp(matom1,matom2,i_pair,jspin))*(conjg(abG1(ll1atom1+1,iG1))*abG2(ll1atom2+1,iG2) &
                                + conjg(abG1(ll1atom1+1+abSizeG1,iG1))*abG2(ll1atom2+1,iG2)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5) &
                                + conjg(abG1(ll1atom1+1,iG1))*abG2(ll1atom2+abSizeG2,iG2)*(usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5) &
                                + conjg(abG1(ll1atom1+1+abSizeG1,iG1))*abG2(ll1atom2+1+abSizeG2,iG2)*(usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5)&
                                *(usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5))*&
                                EXP(cmplx(0.0,-tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),(kpts%bk(:,kptindx)+lapw%gvec(3, iG2,jspin))))
                            ENDDO
                        ENDDO
                        IF(hmat%l_real) THEN
                           hmat%data_r(iG1,iG2)=hmat%data_r(iG1,iG2)+REAL(c_0)
                        ELSE
                           hmat%data_c(iG1,iG2)=hmat%data_c(iG1,iG2)+c_0
                        END IF
                        WRITE(4000,*) 'i_pair,G1,G2,c_0',i_pair,iG1,iG2,c_0
                    ENDDO
                ENDDO
                i_pair=i_pair+1
            ENDDO
        ENDDO
    END SUBROUTINE v_ham
END MODULE m_vham 
