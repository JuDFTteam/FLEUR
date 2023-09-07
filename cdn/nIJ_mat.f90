MODULE m_nIJmat
     
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
                             c_0 = c_0 + we(i) * (conjg(A2)*A1 + conjg(A2)*B1*norm1_W + conjg(B2)*A1*norm2_W + conjg(B2)*B1*norm1_W*norm2_W) * power_factor * exponent
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




        !counter=0
        !DO i_v = 1,atoms%n_v
        !    Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
        !        counter=counter+1
        !    ENDDO
        !ENDDO
        !ALLOCATE(contribution(counter))


        !ALLOCATE(c_0_pair(counter,MAXVAL(lapw%nv),MAXVAL(lapw%nv)))
        !i_pair=1 
        !DO i_v = 1,atoms%n_v  
        !    V_inp=atoms%lda_v(i_v)%V
        !    natom1=atoms%lda_v(i_v)%atomIndex
        !    latom1=atoms%lda_v(i_v)%thisAtomL
        !    ll1atom1=latom1*(latom1+1)
        !    norm1_W = usdus%ddn(latom1,atoms%itype(natom1),jspin)**0.5
        !    CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom1),jspin)
        !    CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom1),natom1,cell,lapw,fjgj,abG1,abSizeG1,.FALSE.)
        !    Do atom2=1,atoms%lda_v(i_v)%numOtherAtoms
        !        natom2=atoms%lda_v(i_v)%otherAtomIndices(atom2)
        !        latom2=atoms%lda_v(i_v)%otherAtomL
        !        ll1atom2=latom2*(latom2+1)
        !        norm2_W = usdus%ddn(latom2,atoms%itype(natom2),jspin)**0.5
        !        power_fac=(cmplx(0, -1)**latom1) *(cmplx(0, 1)**latom2) 
        !        CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,atoms%itype(natom2),jspin)
        !        CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,1,atoms%itype(natom2),natom2,cell,lapw,fjgj,abG2,abSizeG2,.FALSE.)
        !        T5=0.0
        !        T6=0.0
        !        call cpu_time(T5)
        !        DO iG1=1,lapw%nv(jspin)
        !            WRITE(505,*) ' G1, G1 index', iG1, abG1(:,iG1), SIZE(abG1(:,iG1)),abSizeG1
        !            Do iG2=1,lapw%nv(jspin)
        !                c_0=cmplx_0
        !                exponent=EXP(cmplx(0.0,tpi_const)*dot_product(atoms%lda_v(i_v)%atomShifts(:,atom2),(kpts%bk(:,kptindx)+lapw%gvec(:, iG2,jspin))))
        !                Do matom1=-latom1,latom1  
        !                    lm1atom1=ll1atom1+matom1
        !                    Do matom2=-latom2,latom2
        !                        lm1atom2=ll1atom2+matom2
        !                        a1      = abG1(lm1atom1+1,iG1)
        !                        b1      = abG1(lm1atom1+1+abSizeG1/2,iG1)
        !                        a2      = abG2(lm1atom2+1,iG2)
        !                        b2      = abG2(lm1atom2+1+abSizeG2/2,iG2)
        !                        c_0=c_0 - (V_inp)*(conjg(den%nIJ_llp_mmp(matom1,matom2,i_pair,jspin)))* exponent *power_fac * &
        !                        (conjg(a1)*a2 + conjg(b1)*a2*norm1_W + conjg(a1)*b2*norm2_W + conjg(b1)*b2*norm2_W*norm1_W) 
        !                    ENDDO
        !                ENDDO
        !                c_0_pair(i_pair,iG1,iG2)=c_0 
        !            ENDDO
        !        ENDDO
        !        call cpu_time(T6)
        !        WRITE(20,*) 'intenal loop call', T6-T5
        !        WRITE(20,*) 'pair index', i_pair
        !        !! CALL zherk("U","C",lapw%nv(jspin)*lapw%nv(jspin),lapw%nv(jspin),cmplx(1.,0),array(1,lapw%nv(jspin)),1,cmplx(0.,0),lapw%nv(jspin)*lapw%nv(jspin),lapw%nv(jspin))
        !        !! CALL zherk("U","C",lapw%nv(jspin),abSizeG1,cmplx(1.,0),array(1,lapw%nv(jspin)),1,cmplx(0.,0),lapw%nv(jspin)*lapw%nv(jspin),lapw%nv(jspin))
        !        i_pair=i_pair+1
        !    ENDDO
        !ENDDO

        !T7=0.0
        !T8=0.0
        !call cpu_time(T7)
        !DO iG1=1,lapw%nv(jspin)
        !    Do iG2=1,lapw%nv(jspin)
        !        c_pair=cmplx_0
        !        DO indx_pair=1, i_pair-1
        !            c_pair=c_pair+c_0_pair(indx_pair,iG1,iG2)
        !        ENDDO
        !        IF(hmat%l_real) THEN
        !            hmat%data_r(iG1,iG2)=hmat%data_r(iG1,iG2)+REAL(c_pair)
        !        ELSE
        !            hmat%data_c(iG1,iG2)=hmat%data_c(iG1,iG2)+c_pair
        !        END IF
        !    ENDDO
        !ENDDO
        !call cpu_time(T8)
        !WRITE(20,*) 'result loop call', T8-T7



                 !Vectora1= conjg(abG1(ll1atom1-latom1+1:ll1atom1+latom1+1,:))
                 !Vectorb1= conjg(abG1(ll1atom1-latom1+abSizeG1/2+1:ll1atom1+latom1+abSizeG1/2+1,:))
                 !Vectora2= abG2(ll1atom2-latom2+1:ll1atom2+latom2+1,:)
                 !Vectorb2= abG2(ll1atom2-latom2+abSizeG2/2+1:ll1atom2+latom2+abSizeG2/2+1,:)
                 !Matrix_nIJ(:,:,i_pair) = conjg(den%nIJ_llp_mmp(-latom1:latom1,-latom2:latom2,i_pair,jspin))
