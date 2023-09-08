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
        CLASS(t_mat),        INTENT(INOUT)  :: hmat     


        INTEGER i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,iG2,abSizeG1,abSizeG2,counter, atom1index, atom2index
        COMPLEX  power_fac, exponent
        REAL norm1_W, norm2_W, V_inp
        COMPLEX, ALLOCATABLE :: abG1(:,:),abG2(:,:),X1(:,:), X2(:,:), PotMat(:,:)

        ALLOCATE(abG1(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
        ALLOCATE(abG2(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))

        counter=0
        atom1index=0
        atom2index=0
        DO i_v = 1,atoms%n_v
            Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                counter=counter+1
                atom1index= atom1index + 2 * ( 2 * atoms%lda_v(i_v)%thisAtomL +1 )
                atom2index= atom2index + 2 * ( 2 * atoms%lda_v(i_v)%otherAtomL +1 )
            ENDDO
        ENDDO

        ALLOCATE(X1(atom1index,lapw%nv(jspin)))
        ALLOCATE(X2(atom2index,lapw%nv(jspin)))
        ALLOCATE(PotMat(atom1index,atom2index))
        !WRITE(12123,*)'lapw%nv(jspin), MAXVAL(lapw%nv)', lapw%nv(jspin), MAXVAL(lapw%nv)
        i_pair=1 
        atom1index=0
        atom2index=0
        DO i_v = 1,atoms%n_v
            V_inp=atoms%lda_v(i_v)%V
            natom1=atoms%lda_v(i_v)%atomIndex
            latom1=atoms%lda_v(i_v)%thisAtomL
            ll1atom1=latom1*(latom1+1)
            norm1_W = usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5
            CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom1),jspin)
            CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom1),natom1,cell,lapw,fjgj,abG1,abSizeG1,.FALSE.)
            DO atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
                latom2=atoms%lda_v(i_v)%otherAtomL
                ll1atom2=latom2*(latom2+1)
                norm2_W = usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5
                power_fac=(cmplx(0, -1)**latom1) *(cmplx(0, 1)**latom2) 
                PotMat(1+atom1index:atom1index+(2*latom1 +1),1+atom2index:atom2index+(2*latom2 +1)) = conjg(den%nIJ_llp_mmp(-latom1:latom1,-latom2:latom2,i_pair,jspin))
                PotMat(1+atom1index+(2*latom1 +1):atom1index+2*(2*latom1 +1),1+atom2index+(2*latom2 +1):atom2index+2*(2*latom2 +1)) = conjg(den%nIJ_llp_mmp(-latom1:latom1,-latom2:latom2,i_pair,jspin))
                !WRITE(1212,*) 'PotMat',PotMat, i_pair
                CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom2),jspin)
                CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom2),natom2,cell,lapw,fjgj,abG2,abSizeG2,.FALSE.)
                X1( atom1index +1 : atom1index + 2 * latom1 + 1,:) = abG1(ll1atom1-latom1+1:ll1atom1+latom1+1,:)*power_fac
                X1( atom1index + 2*latom1 +1 +1: atom1index + 2*(2*latom1 +1) ,:)= abG1(ll1atom1-latom1+abSizeG1/2+1:ll1atom1+latom1+abSizeG1/2+1,:)*norm1_W*power_fac
                atom1index = atom1index + 2*(2*latom1 +1)
                DO iG2= 1,lapw%nv(jspin) 
                    exponent=EXP(cmplx(0.0,tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),(kpts%bk(:,kptindx)+lapw%gvec(:, iG2,jspin))))
                    X2( atom2index +1 : atom2index + 2 * latom2 + 1,iG2) = abG2(ll1atom2-latom2+1:ll1atom2+latom2+1,iG2)*power_fac*exponent
                    X2( atom2index + 2*latom2 +1 +1: atom2index + 2*(2*latom2 +1) ,iG2)= abG2(ll1atom2-latom2+abSizeG2/2+1:ll1atom2+latom2+abSizeG2/2+1,iG2)*norm2_W*power_fac*exponent
                ENDDO
                atom2index = atom2index + 2*(2*latom2 +1)
                i_pair=i_pair +1
            ENDDO
        ENDDO  
        WRITE(12123,*) 'X2', X2
        WRITE(12123,*) 'potmatr', PotMat
        WRITE(1212,*) 'dimensions',SIZE(X2), SIZE(PotMat), SHAPE(X2), SHAPE(PotMat), RANK(X2), RANK(PotMat)
        WRITE(56563,*) 'sss', MATMUL(PotMat,X2)
        !WRITE(12123,*) 'ssslmllm', MATMUL ( TRANSPOSE(conjg(X1)), MATMUL(PotMat,X2) ) 
        !hmat = hmat + MATMUL(TRANSPOSE(conjg(X1)), MATMUL(PotMat,X2))
        !WRITE(12123,*) 'hmat', SHAPE(MATMUL(TRANSPOSE(conjg(X1)), MATMUL(PotMat,X2)))
        !IF(hmat%l_real) THEN
        !    hmat%data_r=hmat%data_r + MATMUL(TRANSPOSE(conjg(X1)), MATMUL(PotMat,X2))
        !ELSE
        !    hmat%data_c=hmat%data_c + MATMUL(TRANSPOSE(conjg(X1)), MATMUL(PotMat,X2))
        !END IF
        WRITE(4545,*)'hmat'
    END SUBROUTINE v_ham
END MODULE m_vham
