!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP not_used
#else
#define CPP_OMP $OMP
#endif
MODULE m_hsmt_soc_offdiag
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
#ifdef _OPENACC
  SUBROUTINE hsmt_soc_offdiag(n,atoms,cell,fmpi,nococonv,lapw,sym,usdus,td,fjgj,hmat)
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_hsmt_spinor
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_sym  ),INTENT(IN)      :: sym
    TYPE(t_usdus),INTENT(IN)      :: usdus
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_fjgj),INTENT(IN)       :: fjgj
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)!(2,2)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n
    !     ..
    !     ..
    !     .. Local Scalars ..
    INTEGER kii,ki,kj,l,nn,j1,j2,l3
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
    COMPLEX:: chi(2,2,2,2)
    REAL :: plegend(0:2),dplegend(0:2),cross_k(3)
    REAL :: xlegend, dot
    COMPLEX :: cph,fct(2,2),isigma(2,2,3)

    CALL timestart("offdiagonal soc-setup")

    associate(h11=>hmat(1,1)%data_c,h12=>hmat(1,2)%data_c,h21=>hmat(2,1)%data_c,h22=>hmat(2,2)%data_c)
                
    DO l = 0,atoms%lmaxd
       fleg1(l) = REAL(l+l+1)/REAL(l+1)
       fleg2(l) = REAL(l)/REAL(l+1)
       fl2p1(l) = REAL(l+l+1)/fpi_const
    END DO
    !Set up spinors...
    CALL hsmt_spinor_soc(n,nococonv,chi,isigma)
       
    !$acc parallel present(fjgj,fjgj%fj,fjgj%gj,h11,h12,h21,h22)create(cph,dot,fct,xlegend,plegend,dplegend) default(none) copyin(n) &
    !$acc &copyin(chi,isigma,td,td%rsoc,td%rsoc%rso)&
    !$acc &copyin(nococonv,nococonv%alph,nococonv%beta,lapw,lapw%nv,lapw%gvec,lapw%gk,atoms,atoms%firstatom,atoms%neq,atoms%taual,atoms%lmax)&
    !$acc &copyin(fmpi,fleg1,fleg2,fl2p1) 
        
    !$acc loop independent gang private(kii,ki)
    DO  ki =  fmpi%n_rank+1, lapw%nv(1), fmpi%n_size
       kii=(ki-1)/fmpi%n_size+1      
       !$acc loop vector private(cph,dot,fct,xlegend,plegend,dplegend,nn,kj,cross_k)
       DO  kj = 1, ki
        cross_k(1)=lapw%gk(2,ki,1)*lapw%gk(3,kj,1)- lapw%gk(3,ki,1)*lapw%gk(2,kj,1)
        cross_k(2)=lapw%gk(3,ki,1)*lapw%gk(1,kj,1)- lapw%gk(1,ki,1)*lapw%gk(3,kj,1)
        cross_k(3)=lapw%gk(1,ki,1)*lapw%gk(2,kj,1)- lapw%gk(2,ki,1)*lapw%gk(1,kj,1)
           
          !---> set up phase factors
          cph = 0.0          
          DO nn = atoms%firstAtom(n), atoms%firstAtom(n) + atoms%neq(n) - 1
             dot = tpi_const*DOT_PRODUCT(lapw%gvec(:,ki,1)-lapw%gvec(:,kj,1),atoms%taual(:,nn))
             cph = cph + CMPLX(COS(dot),SIN(dot))
          END DO
          !--->       x for legendre polynomials
          xlegend = DOT_PRODUCT(lapw%gk(1:3,kj,1),lapw%gk(1:3,ki,1))
          plegend(0) = 1.0
          dplegend(0) = 0.0
          !--->          update overlap and l-diagonal hamiltonian matrix
          fct=0.0
          DO  l = 1,atoms%lmax(n)
             !--->       legendre polynomials
             l3 = MODULO(l, 3)
             IF (l == 1) THEN
                plegend(1) = xlegend
                dplegend(1) = 1.0
             ELSE
                plegend(l3) = fleg1(l-1)*xlegend*plegend(MODULO(l-1,3)) - fleg2(l-1)*plegend(MODULO(l-2,3))
                dplegend(l3)=REAL(l)*plegend(MODULO(l-1,3))+xlegend*dplegend(MODULO(l-1,3))
             END IF ! l
             DO j2=1,2
                DO j1=1,2
                    fct(j1,j2)  = fct(j1,j2)+cph * dplegend(l3)*fl2p1(l)*(&
                    fjgj%fj(ki,l,j1,1)*fjgj%fj(kj,l,j2,1) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                    fjgj%fj(ki,l,j1,1)*fjgj%gj(kj,l,j2,1) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                    fjgj%gj(ki,l,j1,1)*fjgj%fj(kj,l,j2,1) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                    fjgj%gj(ki,l,j1,1)*fjgj%gj(kj,l,j2,1) *td%rsoc%rso(2,2,n,l,j1,j2)) &
                    * (isigma(j1,j2,1)*cross_k(1)+isigma(j1,j2,2)*cross_k(2)+ isigma(j1,j2,3)*cross_k(3))
                ENDDO
              ENDDO
          ENDDO ! loop over l
          h11(kj,kii)=h11(kj,kii) + chi(1,1,1,1)*fct(1,1)+chi(1,1,1,2)*fct(1,2)+chi(1,1,2,1)*fct(2,1)+chi(1,1,2,2)*fct(2,2)
          h12(kj,kii)=h12(kj,kii) + chi(1,2,1,1)*fct(1,1)+chi(1,2,1,2)*fct(1,2)+chi(1,2,2,1)*fct(2,1)+chi(1,2,2,2)*fct(2,2)
          h21(kj,kii)=h21(kj,kii) + chi(2,1,1,1)*fct(1,1)+chi(2,1,1,2)*fct(1,2)+chi(2,1,2,1)*fct(2,1)+chi(2,1,2,2)*fct(2,2)
          h22(kj,kii)=h22(kj,kii) + chi(2,2,1,1)*fct(1,1)+chi(2,2,1,2)*fct(1,2)+chi(2,2,2,1)*fct(2,1)+chi(2,2,2,2)*fct(2,2)
       ENDDO ! loop over kj
       !$acc end loop
    ENDDO ! loop over ki
    !CPP_OMP END DO
    !$acc end parallel
    !CPP_OMP END PARALLEL
    end associate
   
   
    CALL timestop("offdiagonal soc-setup")

    if (atoms%nlo(n)>0) THEN
      call hsmt_soc_offdiag_LO(n,atoms,cell,fmpi,nococonv,lapw,sym,td,usdus,fjgj,hmat)
    endif  
    RETURN
  END SUBROUTINE hsmt_soc_offdiag
#else
  SUBROUTINE hsmt_soc_offdiag(n,atoms,cell,fmpi,nococonv,lapw,sym,usdus,td,fjgj,hmat)
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_hsmt_spinor
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_sym  ),INTENT(IN)      :: sym
    TYPE(t_usdus),INTENT(IN)      :: usdus
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_fjgj),INTENT(IN)       :: fjgj
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)!(2,2)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n
    !     ..
    !     ..
    !     .. Local Scalars ..
    REAL tnn(3),ski(3), fjkiln,gjkiln
    INTEGER kii,ki,kj,l,nn,j1,j2,ll,l3,kj_off,kj_vec,jv
    INTEGER NVEC_rem  !remainder
    INTEGER, PARAMETER :: NVEC = 128
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd),cross_k(3)
    COMPLEX:: chi(2,2,2,2),isigma(2,2,3)
    REAL, ALLOCATABLE :: plegend(:,:),dplegend(:,:)
    REAL, ALLOCATABLE :: xlegend(:), dot(:)
    COMPLEX, ALLOCATABLE :: cph(:),fct(:),angso(:,:,:)

    CALL timestart("offdiagonal soc-setup")

    !$acc update self(hmat(1,1)%data_c,hmat(2,1)%data_c,hmat(1,2)%data_c,hmat(2,2)%data_c)

    DO l = 0,atoms%lmaxd
       fleg1(l) = REAL(l+l+1)/REAL(l+1)
       fleg2(l) = REAL(l)/REAL(l+1)
       fl2p1(l) = REAL(l+l+1)/fpi_const
    END DO
    !!$acc data copyin(td,td%rsoc,td%rsoc%rso)
    !CPP_OMP PARALLEL DEFAULT(NONE)&
    !CPP_OMP SHARED(n,lapw,atoms,td,fjgj,nococonv,fl2p1,fleg1,fleg2,hmat,fmpi)&
    !CPP_OMP PRIVATE(kii,ki,ski,kj,plegend,dplegend,l,j1,j2,angso,chi)&
    !CPP_OMP PRIVATE(cph,dot,nn,tnn,fct,xlegend,l3,fjkiln,gjkiln,NVEC_rem)&
    !CPP_OMP PRIVATE(kj_off,kj_vec,jv,cross_k,isigma)
    ALLOCATE(cph(NVEC))
    ALLOCATE(xlegend(NVEC))
    ALLOCATE(plegend(NVEC,0:2))
    ALLOCATE(dplegend(NVEC,0:2))
    ALLOCATE(fct(NVEC))
    ALLOCATE(dot(NVEC))
    ALLOCATE(angso(NVEC,2,2))
    !CPP_OMP DO SCHEDULE(DYNAMIC,1)
    DO  ki =  fmpi%n_rank+1, lapw%nv(1), fmpi%n_size
       kii=(ki-1)/fmpi%n_size+1

       DO  kj_off = 1, ki, NVEC
          NVEC_rem = NVEC
          kj_vec = kj_off - 1 + NVEC
          IF (kj_vec > ki) THEN
             kj_vec = ki
             NVEC_rem = ki - kj_off + 1
          ENDIF
          if (NVEC_rem<0 ) exit

          !Set up spinors...
          CALL hsmt_spinor_soc(n,nococonv,chi,isigma)
          DO jv = 1,NVEC_rem
            kj = kj_off - 1 + jv
            cross_k(1)=lapw%gk(2,ki,1)*lapw%gk(3,kj,1)- lapw%gk(3,ki,1)*lapw%gk(2,kj,1)
            cross_k(2)=lapw%gk(3,ki,1)*lapw%gk(1,kj,1)- lapw%gk(1,ki,1)*lapw%gk(3,kj,1)
            cross_k(3)=lapw%gk(1,ki,1)*lapw%gk(2,kj,1)- lapw%gk(2,ki,1)*lapw%gk(1,kj,1)
            DO j1=1,2
              DO j2=1,2
                angso(jv,j1,j2)= (isigma(j1,j2,1)*cross_k(1)+&
                            isigma(j1,j2,2)*cross_k(2)+ isigma(j1,j2,3)*cross_k(3))
              ENDDO
            ENDDO
           ENDDO
       

          !--->             set up phase factors
          cph = 0.0
          ski = lapw%gvec(:,ki,1)
          DO nn = atoms%firstAtom(n), atoms%firstAtom(n) + atoms%neq(n) - 1
             tnn = tpi_const*atoms%taual(:,nn)
             DO jv = 1,NVEC_rem
                kj = kj_off - 1 + jv
                dot(jv) = DOT_PRODUCT(ski(1:3)-lapw%gvec(1:3,kj,1),tnn(1:3))
             END DO
             cph(:NVEC_rem) = cph(:NVEC_rem) + CMPLX(COS(dot(:NVEC_rem)),SIN(dot(:NVEC_rem)))
          END DO

          !--->       x for legendre polynomials
          DO jv = 1,NVEC_rem
             kj = kj_off - 1 + jv
             xlegend(jv) = DOT_PRODUCT(lapw%gk(1:3,kj,1),lapw%gk(1:3,ki,1))
          END DO
          plegend(:NVEC_rem,0) = 1.0
          dplegend(:NVEC_rem,0) = 0.0

          !--->          update overlap and l-diagonal hamiltonian matrix
          !!$acc kernels &
          !!$acc copyin(atoms,atoms%lmax,xlegend,cph,angso)&
          !!$acc create(plegend,dplegend,fct)&
          !!$acc present(fjgj,fjgj%fj,fjgj%gj)&
          !!$acc present(hmat(1,1)%data_c,hmat(2,1)%data_c,hmat(1,2)%data_c,hmat(2,2)%data_c)
          DO  l = 1,atoms%lmax(n)
             !--->       legendre polynomials
             l3 = MODULO(l, 3)
             IF (l == 1) THEN
                plegend(:NVEC_rem,1) = xlegend(:NVEC_rem)
                dplegend(:NVEC_rem,1) = 1.0
             ELSE
                plegend(:NVEC_rem,l3) = fleg1(l-1)*xlegend(:NVEC_rem)*plegend(:NVEC_rem,MODULO(l-1,3)) - fleg2(l-1)*plegend(:NVEC_rem,MODULO(l-2,3))
                dplegend(:NVEC_rem,l3)=REAL(l)*plegend(:NVEC_rem,MODULO(l-1,3))+xlegend(:NVEC_rem)*dplegend(:NVEC_rem,MODULO(l-1,3))
             END IF ! l
             DO j1=1,2
                DO j2=1,2      
                  fct(:NVEC_rem)  =cph(:NVEC_rem) * dplegend(:NVEC_rem,l3)*fl2p1(l)*(&
                  fjgj%fj(ki,l,j1,1)*fjgj%fj(kj_off:kj_vec,l,j2,1) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                  fjgj%fj(ki,l,j1,1)*fjgj%gj(kj_off:kj_vec,l,j2,1) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                  fjgj%gj(ki,l,j1,1)*fjgj%fj(kj_off:kj_vec,l,j2,1) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                  fjgj%gj(ki,l,j1,1)*fjgj%gj(kj_off:kj_vec,l,j2,1) *td%rsoc%rso(2,2,n,l,j1,j2)) &
                  * angso(:NVEC_rem,j1,j2)

                  hmat(1,1)%data_c(kj_off:kj_vec,kii)=hmat(1,1)%data_c(kj_off:kj_vec,kii) + chi(1,1,j1,j2)*fct(:NVEC_rem)
                  hmat(1,2)%data_c(kj_off:kj_vec,kii)=hmat(1,2)%data_c(kj_off:kj_vec,kii) + chi(1,2,j1,j2)*fct(:NVEC_rem)
                  hmat(2,1)%data_c(kj_off:kj_vec,kii)=hmat(2,1)%data_c(kj_off:kj_vec,kii) + chi(2,1,j1,j2)*fct(:NVEC_rem)
                  hmat(2,2)%data_c(kj_off:kj_vec,kii)=hmat(2,2)%data_c(kj_off:kj_vec,kii) + chi(2,2,j1,j2)*fct(:NVEC_rem)
                ENDDO
             ENDDO
          !--->          end loop over l
          ENDDO
          !!$acc end kernels
       ENDDO
    !--->    end loop over ki
    ENDDO
    !CPP_OMP END DO
    !--->       end loop over atom types (ntype)
    DEALLOCATE(xlegend,plegend,dplegend)
    DEALLOCATE(cph)
    !CPP_OMP END PARALLEL
    !!$acc end data
    CALL timestop("offdiagonal soc-setup")

    if (atoms%nlo(n)>0) call hsmt_soc_offdiag_LO(n,atoms,cell,fmpi,nococonv,lapw,sym,td,usdus,fjgj,hmat)
    !$acc update device(hmat(1,1)%data_c,hmat(2,1)%data_c,hmat(1,2)%data_c,hmat(2,2)%data_c)
    RETURN
  END SUBROUTINE hsmt_soc_offdiag


#endif  
  SUBROUTINE hsmt_soc_offdiag_LO(n,atoms,cell,fmpi,nococonv,lapw,sym,td,ud,fjgj,hmat)
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_hsmt_spinor
    USE m_setabc1lo
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_usdus),INTENT(IN)      :: ud
    TYPE(t_fjgj),INTENT(IN)       :: fjgj
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)!(2,2)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n
    !     ..
    !     ..
    !     .. Local Scalars ..
    REAL tnn(3),ski(3)
    INTEGER ki,kj,l,nn,j1,j2,lo,ilo,locol_loc,locol_mat,lorow,ll,nkvec,nkvecp,na,invsfct
    COMPLEX :: fct
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd),cross_k(3)
    COMPLEX:: chi(2,2,2,2),isigma(2,2,3)
    REAL, ALLOCATABLE :: plegend(:,:),dplegend(:,:)
    COMPLEX, ALLOCATABLE :: cph(:)
    REAL                 :: alo1(atoms%nlod,2),blo1(atoms%nlod,2),clo1(atoms%nlod,2)
    INTEGER              :: lo_slot(atoms%nlod),lo_cnt(0:atoms%lmaxd)
    CALL timestart("offdiagonal soc-setup LO")

    DO l = 0,atoms%lmaxd
      fleg1(l) = REAL(l+l+1)/REAL(l+1)
      fleg2(l) = REAL(l)/REAL(l+1)
      fl2p1(l) = REAL(l+l+1)/fpi_const
    END DO
    ALLOCATE(cph(MAXVAL(lapw%nv)))
    ALLOCATE(plegend(MAXVAL(lapw%nv),0:atoms%lmaxd))
    ALLOCATE(dplegend(MAXVAL(lapw%nv),0:atoms%lmaxd))
    plegend=0.0
    plegend(:,0)=1.0
    dplegend(:,0)=0.e0
    dplegend(:,1)=1.e0

    DO j1=1,2
      call setabc1lo(atoms,n,ud,j1, alo1,blo1,clo1)
    ENDDO
    !Normalization taken from hsmt_ab
    alo1=alo1*fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)
    blo1=blo1*fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)
    clo1=clo1*fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)

    !Map each LO to its radial-function slot in rsoc%rso: slot 1=u, 2=udot,
    !3.. = LOs of the same l in the order they appear in atoms%llo (same ordering
    !as in types_radfun%generate_radial_functions).
    lo_cnt = 0
    DO lo = 1,atoms%nlo(n)
       l = atoms%llo(lo,n)
       lo_cnt(l) = lo_cnt(l) + 1
       lo_slot(lo) = 2 + lo_cnt(l)
    ENDDO

    associate(h11=>hmat(1,1)%data_c,h12=>hmat(1,2)%data_c,h21=>hmat(2,1)%data_c,h22=>hmat(2,2)%data_c)

    DO na = atoms%firstAtom(n), atoms%firstAtom(n) + atoms%neq(n) - 1
      IF ((sym%invsat(na) == 0) .OR. (sym%invsat(na) == 1)) THEN
        !--->    if this atom is the first of two atoms related by inversion,
        !--->    the contributions to the overlap matrix of both atoms are added
        !--->    at once. where it is made use of the fact, that the sum of
        !--->    these contributions is twice the real part of the contribution
        !--->    of each atom. note, that in this case there are twice as many
        !--->    (2*(2*l+1)) k-vectors (compare abccoflo and comments there).
        IF (sym%invsat(na) == 0) invsfct = 1
        IF (sym%invsat(na) == 1) invsfct = 2
        !
        DO lo = 1,atoms%nlo(n)
          l = atoms%llo(lo,n)
          if (l==0) cycle !no SOC for s-states
          DO nkvec = 1,invsfct* (2*l+1)
            locol_mat= lapw%nv(1)+lapw%index_lo(lo,na)+nkvec !this is the column of the matrix
            IF (MOD(locol_mat-1,fmpi%n_size) == fmpi%n_rank) THEN !only this MPI rank calculates this column
              locol_loc=(locol_mat-1)/fmpi%n_size+1 !this is the column in local storage
              ki=lapw%kvec(nkvec,lo,na) !this LO is attached to this k+G

              !--->       legendre polynomials
              DO kj = 1,lapw%nv(1)
                plegend(kj,1) = DOT_PRODUCT(lapw%gk(1:3,kj,1),lapw%gk(1:3,ki,1))
              END DO
              DO ll = 1,l - 1
                plegend(:,ll+1) = fleg1(l)*plegend(:,1)*plegend(:,ll) - fleg2(l)*plegend(:,ll-1)
                dplegend(:,ll+1)=REAL(l+1)*plegend(:,ll)+plegend(:,1)*dplegend(:,ll)
              END DO
              !--->             set up phase factors
              cph = 0.0
              ski = lapw%gvec(:,ki,1)
              tnn = tpi_const*atoms%taual(:,na)
              DO kj = 1,lapw%nv(1)
                cph(kj) = cph(kj) +&
                CMPLX(COS(DOT_PRODUCT(ski-lapw%gvec(:,kj,1),tnn)),&
                SIN(DOT_PRODUCT(ski-lapw%gvec(:,kj,1),tnn)))
              END DO
              !Set up spinors...
              CALL hsmt_spinor_soc(n,nococonv,chi,isigma)


              !$acc kernels default(none) &
              !$acc &present(h11,h12,h21,h22)&
              !$acc &present(fjgj,fjgj%fj,fjgj%gj)&
              !$acc &copyin(chi,isigma,td,td%rsoc,td%rsoc%rso)&
              !$acc &copyin(lapw,lapw%gk,lapw%nv,lapw%index_lo,lapw%kvec)&
              !$acc &copyin(alo1,blo1,clo1,cph,dplegend,fl2p1,atoms,atoms%nlo,atoms%llo,lo_slot)&
              !$acc &create(cross_k) 
              DO j1=1,2
                DO j2=1,2
                  !DO j2=j1,j1
                  !---> update l-diagonal hamiltonian matrix with LAPW,LO contribution
                  !$acc loop vector independent private(kj,fct,cross_k)
                  DO kj = 1,lapw%nv(j2)
                    cross_k(1)=lapw%gk(2,ki,1)*lapw%gk(3,kj,1)- lapw%gk(3,ki,1)*lapw%gk(2,kj,1)
                    cross_k(2)=lapw%gk(3,ki,1)*lapw%gk(1,kj,1)- lapw%gk(1,ki,1)*lapw%gk(3,kj,1)
                    cross_k(3)=lapw%gk(1,ki,1)*lapw%gk(2,kj,1)- lapw%gk(2,ki,1)*lapw%gk(1,kj,1)
                    fct  =cph(kj) * dplegend(kj,l)*fl2p1(l)*(&
                    alo1(lo,j1)*fjgj%fj(kj,l,j2,1) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                    alo1(lo,j1)*fjgj%gj(kj,l,j2,1) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                    blo1(lo,j1)*fjgj%fj(kj,l,j2,1) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                    blo1(lo,j1)*fjgj%gj(kj,l,j2,1) *td%rsoc%rso(2,2,n,l,j1,j2)+ &
                    clo1(lo,j1)*fjgj%fj(kj,l,j2,1) *td%rsoc%rso(1,lo_slot(lo),n,l,j1,j2) + &
                    clo1(lo,j1)*fjgj%gj(kj,l,j2,1) *td%rsoc%rso(2,lo_slot(lo),n,l,j1,j2)) &
                    *  (isigma(j1,j2,1)*cross_k(1)+isigma(j1,j2,2)*cross_k(2)+ isigma(j1,j2,3)*cross_k(3))
                    h11(kj,locol_loc)=h11(kj,locol_loc) + chi(1,1,j1,j2)*fct
                    h12(kj,locol_loc)=h12(kj,locol_loc) + chi(1,2,j1,j2)*fct
                    h21(kj,locol_loc)=h21(kj,locol_loc) + chi(2,1,j1,j2)*fct
                    h22(kj,locol_loc)=h22(kj,locol_loc) + chi(2,2,j1,j2)*fct
                  ENDDO
                  !$acc end loop
                  !Update LO-LO part
                  DO ilo=1,atoms%nlo(n)
                    if (l == atoms%llo(ilo,n)) THEN !LO with same L found....
                      DO nkvecp = 1,invsfct* (2*l+1)
                        kj=lapw%kvec(nkvecp,ilo,na) !this LO is attached to this k+G
                        cross_k(1)=lapw%gk(2,ki,1)*lapw%gk(3,kj,1)- lapw%gk(3,ki,1)*lapw%gk(2,kj,1)
                        cross_k(2)=lapw%gk(3,ki,1)*lapw%gk(1,kj,1)- lapw%gk(1,ki,1)*lapw%gk(3,kj,1)
                        cross_k(3)=lapw%gk(1,ki,1)*lapw%gk(2,kj,1)- lapw%gk(2,ki,1)*lapw%gk(1,kj,1)
                        lorow= lapw%nv(1)+lapw%index_lo(ilo,na)+nkvecp !local row
                        if (lorow>locol_mat) cycle
                        fct  =cph(kj) * dplegend(kj,l)*fl2p1(l)*(&
                        alo1(lo,j1)*alo1(ilo,j2) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                        alo1(lo,j1)*blo1(ilo,j2) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                        alo1(lo,j1)*clo1(ilo,j2) *td%rsoc%rso(lo_slot(ilo),1,n,l,j1,j2) + &
                        blo1(lo,j1)*alo1(ilo,j2) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                        blo1(lo,j1)*blo1(ilo,j2) *td%rsoc%rso(2,2,n,l,j1,j2)+ &
                        blo1(lo,j1)*clo1(ilo,j2) *td%rsoc%rso(lo_slot(ilo),2,n,l,j1,j2)+ &
                        clo1(lo,j1)*alo1(ilo,j2) *td%rsoc%rso(1,lo_slot(lo),n,l,j1,j2) + &
                        clo1(lo,j1)*blo1(ilo,j2) *td%rsoc%rso(2,lo_slot(lo),n,l,j1,j2)+ &
                        clo1(lo,j1)*clo1(ilo,j2) *td%rsoc%rso(lo_slot(lo),lo_slot(ilo),n,l,j1,j2)) &
                       *  (isigma(j1,j2,1)*cross_k(1)+isigma(j1,j2,2)*cross_k(2)+ isigma(j1,j2,3)*cross_k(3))
                        h11(lorow,locol_loc)=h11(lorow,locol_loc) + chi(1,1,j1,j2)*fct
                        h12(lorow,locol_loc)=h12(lorow,locol_loc) + chi(1,2,j1,j2)*fct
                        h21(lorow,locol_loc)=h21(lorow,locol_loc) + chi(2,1,j1,j2)*fct
                        h22(lorow,locol_loc)=h22(lorow,locol_loc) + chi(2,2,j1,j2)*fct
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
              !$acc end kernels
            ENDIF !This PE works on LO
          ENDDO!LO
          !--->    end loop over ki
        ENDDO
      ENDIF
    ENDDO
    end associate
    CALL timestop("offdiagonal soc-setup LO")

    RETURN
  END SUBROUTINE hsmt_soc_offdiag_LO
END MODULE m_hsmt_soc_offdiag

#if false
   !Code snipplet useful for debugging only
      !            fct(:NVEC_rem)  =cph(:NVEC_rem) * dplegend(:NVEC_rem,l3)*fl2p1(l)*(&
      !            fjkiln*fjgj%fj(kj_off:kj_vec,l,j2,1) *td%rsoc%rsopp(n,l,j1,j2) ) &
      !            * angso(:NVEC_rem,j1,j2) 
             
      !            BLOCK
      !              use m_anglso
      !              USE m_ylm
      !              INTEGER :: m1,m2,is1,is2,lm1,lm2
      !              COMPLEX :: soangl(0:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,-atoms%lmaxd:atoms%lmaxd,2),angso2
      !              COMPLEX :: ylm1( (atoms%lmaxd+1)**2 ), ylm2( (atoms%lmaxd+1)**2 )
      !              INTEGER :: ispjsp(2)
      !              if (kj_off/=kj_vec) call judft_error("DEBUG Problem")
      !              DATA ispjsp/1,-1/
      !                   CALL ylm4(l,lapw%gk(:,ki,1),ylm1)
      !                   CALL ylm4(l,lapw%gk(:,kj,1),ylm2)
      !                   angso2=0.0
      !                   is1=ispjsp(j1)
      !                   is2=ispjsp(j2)
      !                   DO m1=-l,l
      !                      lm1=l*(l+1)+m1+1
      !                      DO m2=-l,l
      !                        lm2=l*(l+1)+m2+1
      !                        angso2=angso2+ylm1(lm1)*conjg(ylm2(lm2))* &
      !                           anglso(nococonv%beta(n),nococonv%alph(n),l,m1,is1,l,m2,is2)
      !                     ENDDO
      !                  ENDDO
      !                  fct(1)=angso2*fjgj%fj(ki,l,j1,1)*fjgj%fj(kj_off,l,j2,1) *td%rsoc%rsopp(n,l,j1,j2)
      !            END BLOCK     
#endif
