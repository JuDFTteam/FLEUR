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

  !Development switch: set to .TRUE. to have hsmt_soc_offdiag_check verify the
  !closed-form SOC angular factor against an explicit spherical-harmonic reference
  !(see that routine for what is tested). Always compiled, never on in production.
  LOGICAL, PARAMETER :: l_checkSOCangular = .FALSE.

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
                    fjgj%fj(kj,l,j1,1)*fjgj%fj(ki,l,j2,1) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                    fjgj%gj(kj,l,j1,1)*fjgj%fj(ki,l,j2,1) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                    fjgj%fj(kj,l,j1,1)*fjgj%gj(ki,l,j2,1) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                    fjgj%gj(kj,l,j1,1)*fjgj%gj(ki,l,j2,1) *td%rsoc%rso(2,2,n,l,j1,j2)) &
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

    !Development check, off by default; see l_checkSOCangular at the top of this module.
    IF (l_checkSOCangular) CALL hsmt_soc_offdiag_check(n,atoms,fmpi,nococonv,lapw)

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
                  fjgj%fj(kj_off:kj_vec,l,j1,1)*fjgj%fj(ki,l,j2,1) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                  fjgj%gj(kj_off:kj_vec,l,j1,1)*fjgj%fj(ki,l,j2,1) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                  fjgj%fj(kj_off:kj_vec,l,j1,1)*fjgj%gj(ki,l,j2,1) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                  fjgj%gj(kj_off:kj_vec,l,j1,1)*fjgj%gj(ki,l,j2,1) *td%rsoc%rso(2,2,n,l,j1,j2)) &
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
                plegend(:,ll+1) = fleg1(ll)*plegend(:,1)*plegend(:,ll) - fleg2(ll)*plegend(:,ll-1)
                dplegend(:,ll+1)=REAL(ll+1)*plegend(:,ll)+plegend(:,1)*dplegend(:,ll)
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
                    alo1(lo,j2)*fjgj%fj(kj,l,j1,1) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                    alo1(lo,j2)*fjgj%gj(kj,l,j1,1) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                    blo1(lo,j2)*fjgj%fj(kj,l,j1,1) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                    blo1(lo,j2)*fjgj%gj(kj,l,j1,1) *td%rsoc%rso(2,2,n,l,j1,j2)+ &
                    clo1(lo,j2)*fjgj%fj(kj,l,j1,1) *td%rsoc%rso(1,lo_slot(lo),n,l,j1,j2) + &
                    clo1(lo,j2)*fjgj%gj(kj,l,j1,1) *td%rsoc%rso(2,lo_slot(lo),n,l,j1,j2)) &
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
                        alo1(lo,j2)*alo1(ilo,j1) *td%rsoc%rso(1,1,n,l,j1,j2) + &
                        alo1(lo,j2)*blo1(ilo,j1) *td%rsoc%rso(2,1,n,l,j1,j2) + &
                        alo1(lo,j2)*clo1(ilo,j1) *td%rsoc%rso(lo_slot(ilo),1,n,l,j1,j2) + &
                        blo1(lo,j2)*alo1(ilo,j1) *td%rsoc%rso(1,2,n,l,j1,j2) + &
                        blo1(lo,j2)*blo1(ilo,j1) *td%rsoc%rso(2,2,n,l,j1,j2)+ &
                        blo1(lo,j2)*clo1(ilo,j1) *td%rsoc%rso(lo_slot(ilo),2,n,l,j1,j2)+ &
                        clo1(lo,j2)*alo1(ilo,j1) *td%rsoc%rso(1,lo_slot(lo),n,l,j1,j2) + &
                        clo1(lo,j2)*blo1(ilo,j1) *td%rsoc%rso(2,lo_slot(lo),n,l,j1,j2)+ &
                        clo1(lo,j2)*clo1(ilo,j1) *td%rsoc%rso(lo_slot(ilo),lo_slot(lo),n,l,j1,j2)) &
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
  SUBROUTINE hsmt_soc_offdiag_check(n,atoms,fmpi,nococonv,lapw)
    !! Development check for the closed-form SOC angular factor. Always compiled but
    !! only called when l_checkSOCangular (top of this module) is set to .TRUE.; it is
    !! wired into the CPU branch of hsmt_soc_offdiag only.
    !!
    !! hsmt_soc_offdiag evaluates the angular part of the first-variation SOC matrix
    !! element in closed form. With k_row = gk(kj), k_col = gk(ki), x = k_row.k_col and
    !! cross_k = k_col x k_row it uses
    !!
    !!   A_code(j1,j2) = fl2p1(l) * P_l'(x) * sum_a isigma(j1,j2,a) cross_k(a)
    !!
    !! where isigma(j1,j2,a) = <chi_j1|i*sigma_a|chi_j2> in the local spin basis. The
    !! exact expression it stands for is
    !!
    !!   A(j1,j2) = sum_{mr,mc} Y_lmr(k_row) Y*_lmc(k_col) <l mr j1|L.sigma|l mc j2>
    !!            = sum_a V_a(l) <chi_j1|sigma_a|chi_j2>,
    !!   V_a(l)   = sum_{mr,mc} Y_lmr(k_row) Y*_lmc(k_col) <l mr|L_a|l mc>.
    !!
    !! Two independent checks are run:
    !!
    !! 1. "Ylm": the addition theorem  V_a(l) = fl2p1(l) * P_l'(x) * i * cross_k(a),
    !!    evaluated from ylm4 and explicit L_a matrices. Spin-free and valid at any
    !!    spin-frame angle, so it covers P_l' for every l. P_l' is built here with the
    !!    array recursion of hsmt_soc_offdiag_LO, i.e. this is the check that pins the
    !!    LO path (the LAPW-LAPW block uses the equivalent rolling three-slot form).
    !!
    !! 2. "anglso": the full per-spin-block element against ylm4 + anglso, which
    !!    returns <l m1 is1|L.sigma|l m2 is2> directly. This is the only check that
    !!    tests the (m,spin) pairing, i.e. that j1 is the row and j2 the column. It is
    !!    reported together with the row/column-exchanged pairing that preceded this
    !!    convention; the exchanged one must NOT vanish for j1/=j2.
    !!    CAUTION: anglso's default l_standard_euler_angles=.FALSE. branch flips the
    !!    spin x and y axes (see the comment in anglso.f90), so it agrees with the
    !!    standard rotation only for beta=alpha=0, and there only up to a sign in the
    !!    spin-off-diagonal blocks (compensated below). For a rotated frame it is a
    !!    different operator, not a reference, and this check is therefore skipped.
    USE m_constants, ONLY : fpi_const,oUnit
    USE m_types
    USE m_hsmt_spinor
    USE m_anglso
    USE m_ylm
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    INTEGER,INTENT(IN)            :: n

    INTEGER, PARAMETER :: NSAMPLE = 12          !k-vectors sampled per side
    INTEGER, PARAMETER :: ispjsp(2) = [1,-1]
    REAL,    PARAMETER :: TOL = 1.0e-10

    INTEGER :: l,ll,j1,j2,mr,mc,ki,kj,ika,ikb,stride,lmr,lmc,a
    REAL    :: x,dpl,cross_k(3),sgn,lplus,lminus
    REAL    :: fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
    REAL    :: pl(0:atoms%lmaxd),dpleg(0:atoms%lmaxd)
    REAL    :: dev_ylm,scale_ylm,worst,sample_scale,cnorm
    REAL    :: dev(2,2,2),scale_s(2,2)
    LOGICAL :: l_unrot
    COMPLEX :: chi(2,2,2,2),isigma(2,2,3)
    COMPLEX :: a_code,a_ref,a_swap,lel(3),vvec(3),yy
    COMPLEX :: ylm_row((atoms%lmaxd+1)**2),ylm_col((atoms%lmaxd+1)**2)

    IF (fmpi%irank/=0) RETURN
    IF (lapw%nv(1)<2) RETURN

    DO l = 0,atoms%lmaxd
       fleg1(l) = REAL(l+l+1)/REAL(l+1)
       fleg2(l) = REAL(l)/REAL(l+1)
       fl2p1(l) = REAL(l+l+1)/fpi_const
    END DO
    CALL hsmt_spinor_soc(n,nococonv,chi,isigma)

    l_unrot = ABS(nococonv%beta(n))<1.0e-8 .AND. ABS(nococonv%alph(n))<1.0e-8

    WRITE (oUnit,'(/,a,i0,a,f9.5,a,f9.5)') &
         " SOC-angular-check: atom type ",n, &
         "   beta=",nococonv%beta(n),"  alpha=",nococonv%alph(n)
    IF (l_unrot) THEN
       WRITE (oUnit,'(a)') "   columns: |Ylm| = addition-theorem deviation (spin-free);"
       WRITE (oUnit,'(a)') "            |anglso| = full element vs ylm4+anglso, per spin block;"
       WRITE (oUnit,'(a)') "            |swapped| = same with the row/column (m,spin) pairing"
       WRITE (oUnit,'(a)') "            exchanged, which must NOT vanish for j1/=j2."
       WRITE (oUnit,'(a)') "   all deviations relative;  l  j1  j2      |Ylm|    |anglso|   |swapped|"
    ELSE
       WRITE (oUnit,'(a)') "   rotated spin frame: the anglso check is skipped (its default"
       WRITE (oUnit,'(a)') "   l_standard_euler_angles=.FALSE. branch flips the spin x/y axes"
       WRITE (oUnit,'(a)') "   and is not a reference here). Only the spin-free Ylm check runs."
       WRITE (oUnit,'(a)') "   all deviations relative;  l               |Ylm|"
    END IF

    worst = 0.0
    stride = MAX(1,lapw%nv(1)/NSAMPLE)
    DO l = 1,atoms%lmax(n)
       dev = 0.0 ; scale_s = 0.0 ; dev_ylm = 0.0 ; scale_ylm = 0.0
       DO ika = 1,lapw%nv(1),stride
          ki = ika                       !column
          !gk is (k+G)/|k+G|; for k+G=0 (Gamma with G=0) it degenerates to the null
          !vector and carries no direction, so such a basis function is skipped here.
          IF (DOT_PRODUCT(lapw%gk(:,ki,1),lapw%gk(:,ki,1))<0.5) CYCLE
          CALL ylm4(l,lapw%gk(:,ki,1),ylm_col)
          DO ikb = 1,lapw%nv(1),stride
             kj = ikb                    !row
             IF (kj==ki) CYCLE           !cross product vanishes, nothing to test
             IF (DOT_PRODUCT(lapw%gk(:,kj,1),lapw%gk(:,kj,1))<0.5) CYCLE
             CALL ylm4(l,lapw%gk(:,kj,1),ylm_row)

             !--- P_l'(x) with the array recursion of hsmt_soc_offdiag_LO
             x = DOT_PRODUCT(lapw%gk(1:3,kj,1),lapw%gk(1:3,ki,1))
             pl = 0.0 ; dpleg = 0.0
             pl(0) = 1.0 ; dpleg(0) = 0.0
             pl(1) = x   ; dpleg(1) = 1.0
             DO ll = 1,l-1
                pl(ll+1)    = fleg1(ll)*x*pl(ll) - fleg2(ll)*pl(ll-1)
                dpleg(ll+1) = REAL(ll+1)*pl(ll) + x*dpleg(ll)
             END DO
             dpl = dpleg(l)

             cross_k(1)=lapw%gk(2,ki,1)*lapw%gk(3,kj,1)-lapw%gk(3,ki,1)*lapw%gk(2,kj,1)
             cross_k(2)=lapw%gk(3,ki,1)*lapw%gk(1,kj,1)-lapw%gk(1,ki,1)*lapw%gk(3,kj,1)
             cross_k(3)=lapw%gk(1,ki,1)*lapw%gk(2,kj,1)-lapw%gk(2,ki,1)*lapw%gk(1,kj,1)

             !--- check 1: V_a(l) = fl2p1 * P_l' * i * cross_k(a)
             vvec = CMPLX(0.0,0.0)
             DO mr = -l,l
                lmr = l*(l+1)+mr+1
                DO mc = -l,l
                   lmc = l*(l+1)+mc+1
                   yy = ylm_row(lmr)*CONJG(ylm_col(lmc))
                   !<l mr|L_a|l mc>: L+ raises mc, L- lowers it
                   lplus  = 0.0 ; lminus = 0.0
                   IF (mr==mc+1) lplus  = SQRT(REAL((l-mc)*(l+mc+1)))
                   IF (mr==mc-1) lminus = SQRT(REAL((l+mc)*(l-mc+1)))
                   lel(1) = CMPLX(0.5*(lplus+lminus),0.0)              !L_x
                   lel(2) = CMPLX(0.0,-0.5*(lplus-lminus))             !L_y = (L+ - L-)/(2i)
                   lel(3) = CMPLX(0.0,0.0)
                   IF (mr==mc) lel(3) = CMPLX(REAL(mc),0.0)            !L_z
                   vvec(:) = vvec(:) + yy*lel(:)
                END DO
             END DO
             !relative per sample, so a near-degenerate k-pair (|cross_k|->0, where both
             !sides vanish) cannot dominate a global maximum
             sample_scale = 0.0
             DO a = 1,3
                sample_scale = MAX(sample_scale,ABS(vvec(a)))
             END DO
             cnorm = SQRT(DOT_PRODUCT(cross_k,cross_k))
             IF (sample_scale>1.0e-8 .AND. cnorm>1.0e-8) THEN
                DO a = 1,3
                   dev_ylm = MAX(dev_ylm, &
                        ABS(vvec(a)-fl2p1(l)*dpl*CMPLX(0.0,1.0)*cross_k(a))/sample_scale)
                END DO
                scale_ylm = 1.0
             END IF

             !--- check 2: full element against ylm4 + anglso (unrotated frames only)
             IF (.NOT.l_unrot) CYCLE
             DO j1 = 1,2      !row (bra) local spin
                DO j2 = 1,2   !column (ket) local spin
                   a_code = fl2p1(l)*dpl*( isigma(j1,j2,1)*cross_k(1) &
                                          +isigma(j1,j2,2)*cross_k(2) &
                                          +isigma(j1,j2,3)*cross_k(3) )
                   sgn = MERGE(1.0,-1.0,j1==j2)   !anglso convention, see header
                   a_ref  = CMPLX(0.0,0.0)
                   a_swap = CMPLX(0.0,0.0)
                   DO mr = -l,l
                      lmr = l*(l+1)+mr+1
                      DO mc = -l,l
                         lmc = l*(l+1)+mc+1
                         yy = ylm_row(lmr)*CONJG(ylm_col(lmc))
                         a_ref  = a_ref  + yy*anglso(nococonv%beta(n),nococonv%alph(n), &
                                                     l,mr,ispjsp(j1),l,mc,ispjsp(j2))
                         a_swap = a_swap + yy*anglso(nococonv%beta(n),nococonv%alph(n), &
                                                     l,mr,ispjsp(j2),l,mc,ispjsp(j1))
                      END DO
                   END DO
                   dev(1,j1,j2) = MAX(dev(1,j1,j2),ABS(a_code-sgn*a_ref))
                   dev(2,j1,j2) = MAX(dev(2,j1,j2),ABS(a_code-sgn*a_swap))
                   scale_s(j1,j2) = MAX(scale_s(j1,j2),ABS(a_ref))
                END DO
             END DO
          END DO
       END DO

       IF (scale_ylm>1.0e-12) worst = MAX(worst,dev_ylm)
       IF (l_unrot) THEN
          DO j1 = 1,2
             DO j2 = 1,2
                IF (scale_s(j1,j2)<1.0e-12) CYCLE
                worst = MAX(worst,dev(1,j1,j2)/scale_s(j1,j2))
                WRITE (oUnit,'(3i4,3f12.8)') l,j1,j2,dev_ylm, &
                     dev(1,j1,j2)/scale_s(j1,j2),dev(2,j1,j2)/scale_s(j1,j2)
             END DO
          END DO
       ELSE
          WRITE (oUnit,'(i4,12x,f12.8)') l,dev_ylm
       END IF
    END DO

    IF (worst<TOL) THEN
       WRITE (oUnit,'(a,e12.4)') " SOC-angular-check: PASS, worst relative deviation ",worst
    ELSE
       WRITE (oUnit,'(a,e12.4)') " SOC-angular-check: FAIL, worst relative deviation ",worst
    END IF

  END SUBROUTINE hsmt_soc_offdiag_check
END MODULE m_hsmt_soc_offdiag
