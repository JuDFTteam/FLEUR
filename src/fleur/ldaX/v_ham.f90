MODULE m_vham

    !----------------------------------------------------------------------------------- !
    !                                                                                    !
    ! MODULE: m_v_ham                                                                    !
    !                                                                                    !
    ! Author : Wejdan Beida                                                              !
    !                                                                                    !
    ! Description: This module calculates the intersite potential matrix elemnts         !
    !                                                                                    !
    !     s                 s         ---   ---   IJ   JI     s     I     J     s        !
    !<Phsi  |V_POT_MAT| Phis  > = -   >     >    V    n   <Psi  |phi > <phi |Psi  >      !
    !     k                 k         ---   ---               k     L      Lp   k        !
    !                                 I,J,s  L,Lp                                        !
    !                                                                                    !
    ! where L={l,m},Lp={lp,mp}                                                           !
    !                                                                                    !
    !------------------------------------------------------------------------------------! 


    CONTAINS

    SUBROUTINE v_ham(input,usdus,atoms,kpts,cell,lapw,sym,noco,fmpi,nococonv,fjgj,den,jspin,kptindx,hmat)

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
        TYPE(t_mpi),         INTENT(IN)     :: fmpi
        TYPE(t_nococonv),    INTENT(IN)     :: nococonv
        TYPE(t_fjgj),        INTENT(INOUT)  :: fjgj
        TYPE(t_potden),      INTENT(IN)     :: den
        INTEGER,             INTENT(IN)     :: jspin,kptindx
        CLASS(t_mat),        INTENT(INOUT)  :: hmat     




            INTEGER i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2,iG1,iG2,abSizeG1,abSizeG2
            COMPLEX c_0, a1, b1, a2, b2, power_fac, exponent
            REAL norm1_W, norm2_W, V_inp
            COMPLEX, ALLOCATABLE :: abG1(:,:),abG2(:,:),temp_nIJ(:,:), c_pairs(:,:)
            ALLOCATE(abG1(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
            ALLOCATE(abG2(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
            ALLOCATE(temp_nIJ(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
            ALLOCATE(c_pairs(MAXVAL(lapw%nv),MAXVAL(lapw%nv)))

            !c_pairs=cmplx_0
            abG1=cmplx_0
            abG2=cmplx_0
            temp_nIJ=cmplx_0
            c_pairs=cmplx_0


            i_pair=1 
            DO i_v = 1,atoms%n_v  
                V_inp=atoms%lda_v(i_v)%V / hartree_to_ev_const
                natom1=atoms%lda_v(i_v)%atomIndex
                latom1=atoms%lda_v(i_v)%thisAtomL
                ll1atom1=latom1*(latom1+1)
                norm1_W = usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5
                CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom1),jspin)
                CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom1),natom1,cell,lapw,fjgj,abG1,abSizeG1,.FALSE.)
                Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                    natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
                    latom2=atoms%lda_v(i_v)%otherAtomL
                    ll1atom2=latom2*(latom2+1)
                    norm2_W = usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5
                    power_fac=(cmplx(0, -1)**latom1) *(cmplx(0, 1)**latom2) 
                    CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom2),jspin)
                    CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom2),natom2,cell,lapw,fjgj,abG2,abSizeG2,.FALSE.)
                    DO iG2=1,lapw%nv(jspin)
                        exponent=EXP(cmplx(0.0,tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),(kpts%bk(:,kptindx)+lapw%gvec(:, iG2,jspin))))
                        temp_nIJ(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)=TRANSPOSE(conjg(den%nIJ_llp_mmp(:,:,i_pair,jspin)))
                        Do iG1=1,lapw%nv(jspin)
                            c_0=cmplx_0
                            Do matom1=-latom1,latom1  
                                lm1atom1=ll1atom1+matom1
                                a1      = abG1(lm1atom1+1,iG1)
                                b1      = abG1(lm1atom1+1+abSizeG1/2,iG1)
                                Do matom2=-latom2,latom2
                                    lm1atom2=ll1atom2+matom2
                                    a2      = abG2(lm1atom2+1,iG2)
                                    b2      = abG2(lm1atom2+1+abSizeG2/2,iG2)
                                    c_0     = c_0 - (V_inp)*temp_nIJ(matom1,matom2)* exponent * power_fac * &
                                                 (conjg(a1)*a2 + conjg(b1)*a2*norm1_W + conjg(a1)*b2*norm2_W + conjg(b1)*b2*norm2_W*norm1_W)
                                ENDDO
                            ENDDO
                            c_pairs(iG1,iG2)=c_0+c_pairs(iG1,iG2) 
                        ENDDO
                    ENDDO
                    i_pair=i_pair+1
                ENDDO
            ENDDO


            DO iG1=1,lapw%nv(jspin)
                Do iG2=1,lapw%nv(jspin)
                    IF(hmat%l_real) THEN
                        hmat%data_r(iG1,iG2)=hmat%data_r(iG1,iG2)+REAL(c_pairs(iG1,iG2))
                    ELSE
                        hmat%data_c(iG1,iG2)=hmat%data_c(iG1,iG2)+c_pairs(iG1,iG2)
                    END IF
                ENDDO
            ENDDO
    END SUBROUTINE v_ham
END MODULE m_vham
