MODULE m_vham

    !----------------------------------------------------------------------------------- !
    !                                                                                    !
    ! MODULE: m_v_ham                                                                    !
    !                                                                                    !
    ! Author : Wejdan Beida                                                              !
    !                                                                                    !
    ! Description: This module calculates the intersite occupation matrix                !
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




            INTEGER i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,matom1,matom2,lm1atom1,lm1atom2,iG1,iG2,abSizeG1,abSizeG2,indx_pair,counter, atom1index, atom2index
            COMPLEX c_0,c_pair, a1, b1, a2, b2, power_fac, exponent
            REAL norm1_W, norm2_W, V_inp
            COMPLEX, ALLOCATABLE :: abG1(:,:),abG2(:,:)
            COMPLEX, ALLOCATABLE :: c_0_pair(:,:,:)  
            ALLOCATE(abG1(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
            ALLOCATE(abG2(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))

            counter=0
            DO i_v = 1,atoms%n_v
                Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
                    counter=counter+1
                ENDDO
            ENDDO
            ALLOCATE(c_0_pair(counter,MAXVAL(lapw%nv),MAXVAL(lapw%nv)))

            DO iG1=1,lapw%nv(jspin)
                DO iG2=1,lapw%nv(jspin)
                    DO indx_pair=1, counter-1
                        c_0_pair(indx_pair,iG1,iG2)=cmplx_0
                    ENDDO
                ENDDO
            ENDDO

            i_pair=1 
            DO i_v = 1,atoms%n_v  
                V_inp=atoms%lda_v(i_v)%V
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
                    DO iG1=1,lapw%nv(jspin)
                        Do iG2=1,lapw%nv(jspin)
                            c_0=cmplx_0
                            exponent=EXP(cmplx(0.0,tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),(kpts%bk(:,kptindx)+lapw%gvec(:, iG2,jspin))))
                            Do matom1=-latom1,latom1  
                                lm1atom1=ll1atom1+matom1
                                Do matom2=-latom2,latom2
                                    lm1atom2=ll1atom2+matom2
                                    a1      = abG1(lm1atom1+1,iG1)
                                    b1      = abG1(lm1atom1+1+abSizeG1/2,iG1)
                                    a2      = abG2(lm1atom2+1,iG2)
                                    b2      = abG2(lm1atom2+1+abSizeG2/2,iG2)
                                    c_0     = c_0 - (V_inp)*(conjg(den%nIJ_llp_mmp(matom1,matom2,i_pair,jspin)))* exponent *power_fac * &
                                              (conjg(a1)*a2 + conjg(b1)*a2*norm1_W + conjg(a1)*b2*norm2_W + conjg(b1)*b2*norm2_W*norm1_W) 
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
                    DO indx_pair=1, i_pair-1
                        c_pair=c_pair+c_0_pair(indx_pair,iG1,iG2)
                    ENDDO
                    IF(hmat%l_real) THEN
                        hmat%data_r(iG1,iG2)=hmat%data_r(iG1,iG2)+REAL(c_pair)
                    ELSE
                        hmat%data_c(iG1,iG2)=hmat%data_c(iG1,iG2)+c_pair
                    END IF
                ENDDO
            ENDDO


            !INTEGER i_v,i_pair,natom1,latom1,ll1atom1,atom2,natom2,latom2,ll1atom2,iG2,abSizeG1,abSizeG2,counter, atom1index, atom2index, i,j,k
            !COMPLEX  power_fac, exponent
            !REAL norm1_W, norm2_W, V_inp
            !COMPLEX, ALLOCATABLE :: abG1(:,:),abG2(:,:),X1(:,:), X2(:,:), PotMat(:,:)
            !ALLOCATE(abG1(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
            !ALLOCATE(abG2(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
            !atom1index=0
            !atom2index=0
            !counter=0
            !DO i_v = 1,atoms%n_v
            !    Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
            !        counter=counter+1
            !        atom1index= atom1index + 2 * ( 2 * atoms%lda_v(i_v)%thisAtomL +1 )
            !        atom2index= atom2index + 2 * ( 2 * atoms%lda_v(i_v)%otherAtomL +1 )
            !    ENDDO
            !ENDDO
            !ALLOCATE(X1(atom1index,lapw%nv(jspin)))
            !ALLOCATE(X2(atom2index,lapw%nv(jspin)))
            !ALLOCATE(PotMat(atom1index,atom2index))

            !DO i=1,atom1index
            !    DO j=1,atom2index
            !        PotMat(i,j)= cmplx_0
            !        DO k =1, MAXVAL(lapw%nv)
            !            X1(i,k)=cmplx_0
            !            X2(j,k)=cmplx_0
            !        ENDDO
            !    ENDDO
            !ENDDO

            !i_pair=1 
            !atom1index=0
            !atom2index=0
            !counter=0
            !DO i_v = 1,atoms%n_v
            !    V_inp=atoms%lda_v(i_v)%V
            !    natom1=atoms%lda_v(i_v)%atomIndex
            !    latom1=atoms%lda_v(i_v)%thisAtomL
            !    ll1atom1=latom1*(latom1+1)
            !    norm1_W = usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5
            !    CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom1),jspin)
            !    CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom1),natom1,cell,lapw,fjgj,abG1,abSizeG1,.FALSE.)
            !    DO atom2=1,atoms%lda_v(i_v)%numOtherAtoms
            !        natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
            !        latom2=atoms%lda_v(i_v)%otherAtomL
            !        ll1atom2=latom2*(latom2+1)
            !        norm2_W = usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5
            !        power_fac=(cmplx(0, -1)**latom1) *(cmplx(0, 1)**latom2)
            !        X1( atom1index +1 : atom1index + 2 * latom1 + 1,:) = conjg(abG1(ll1atom1-latom1+1:ll1atom1+latom1+1,:))*power_fac*V_inp
            !        X1( atom1index + 2*latom1 +1 +1: atom1index + 2*(2*latom1 +1) ,:)= conjg(abG1(ll1atom1-latom1+abSizeG1/2+1:ll1atom1+latom1+abSizeG1/2+1,:))*norm1_W*power_fac*V_inp
            !        CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom2),jspin)
            !        CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom2),natom2,cell,lapw,fjgj,abG2,abSizeG2,.FALSE.)
            !        DO iG2= 1,lapw%nv(jspin)
            !            exponent=EXP(cmplx(0.0,tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),(kpts%bk(:,kptindx)+lapw%gvec(:, iG2,jspin))))
            !            X2( atom2index +1 : atom2index + 2 * latom2 + 1,:) = abG2(ll1atom2-latom2+1:ll1atom2+latom2+1,:)*exponent
            !            X2( atom2index + 2*latom2 +1 +1: atom2index + 2*(2*latom2 +1) ,:)= abG2(ll1atom2-latom2+abSizeG2/2+1:ll1atom2+latom2+abSizeG2/2+1,:)*exponent*norm2_W
            !        ENDDO
            !        PotMat(1+atom1index:atom1index+(2*latom1 +1),1+atom2index:atom2index+(2*latom2 +1)) = conjg(den%nIJ_llp_mmp(-latom1:latom1,-latom2:latom2,i_pair,jspin))
            !        PotMat(1+atom1index+(2*latom1 +1):atom1index+2*(2*latom1 +1),1+atom2index+(2*latom2 +1):atom2index+2*(2*latom2 +1)) = conjg(den%nIJ_llp_mmp(-latom1:latom1,-latom2:latom2,i_pair,jspin))
            !        atom2index = atom2index + 2*(2*latom2 +1)
            !        atom1index = atom1index + 2*(2*latom1 +1)
            !        i_pair=i_pair +1
            !    ENDDO
            !ENDDO
            !IF(hmat%l_real) THEN
            !    hmat%data_r=hmat%data_r - REAL(MATMUL(TRANSPOSE(X1), MATMUL(TRANSPOSE(PotMat),X2)))
            !ELSE
            !    hmat%data_c=hmat%data_c - MATMUL(TRANSPOSE(X1), MATMUL(TRANSPOSE(PotMat),X2))
            !END IF 

            ! For debugging (END)
            !i_pair=1
            !WRITE(1000+kptindx,*) '=========='
            !DO i_v = 1,atoms%n_v
            !   natom1=atoms%lda_v(i_v)%atomIndex
            !   latom1=atoms%lda_v(i_v)%thisAtomL
            !   Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
            !      natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
            !      latom2=atoms%lda_v(i_v)%otherAtomL
            !      DO matom1=-latom1,latom1
            !         DO matom2=-latom2,latom2
            !            WRITE(1000+kptindx,'(4i7,2f15.8)') matom2,matom1,i_Pair,jspin,den%nIJ_llp_mmp(matom1,matom2,i_pair,jspin)
            !         END DO
            !      END DO
            !      i_pair=i_pair+1
            !   END DO
            !END DO
            ! For debugging (END)

            !        T5=0.0
            !        T6=0.0
            !        call cpu_time(T5)
            !        call cpu_time(T6)
            !         WRITE(20,*) 'intenal loop call', T6-T5
    END SUBROUTINE v_ham
END MODULE m_vham
