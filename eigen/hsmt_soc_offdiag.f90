!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_soc_offdiag
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_soc_offdiag(n,atoms,cell,mpi,nococonv,lapw,sym,usdus,td,fj,gj,hmat)
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_hsmt_spinor
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_sym  ),INTENT(IN)      :: sym
    TYPE(t_usdus),INTENT(IN)      :: usdus
    TYPE(t_tlmplm),INTENT(IN)     :: td
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)!(2,2)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: fj(:,0:,:),gj(:,0:,:)
    !     ..
    !     .. Local Scalars ..
    REAL tnn(3),ski(3)
    INTEGER kii,ki,kj,l,nn,j1,j2,ll
    COMPLEX :: fct
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
    COMPLEX:: chi(2,2,2,2),angso(lapw%nv(1),2,2)
    REAL, ALLOCATABLE :: plegend(:,:),dplegend(:,:)
    COMPLEX, ALLOCATABLE :: cph(:)

    CALL timestart("offdiagonal soc-setup")

    DO l = 0,atoms%lmaxd
       fleg1(l) = REAL(l+l+1)/REAL(l+1)
       fleg2(l) = REAL(l)/REAL(l+1)
       fl2p1(l) = REAL(l+l+1)/fpi_const
    END DO

    !$OMP PARALLEL DEFAULT(NONE)&
    !$OMP SHARED(n,lapw,atoms,td,fj,gj,nococonv,fl2p1,fleg1,fleg2,hmat,mpi)&
    !$OMP PRIVATE(kii,ki,ski,kj,plegend,dplegend,l,j1,j2,angso,chi)&
    !$OMP PRIVATE(cph,nn,tnn,fct)
    ALLOCATE(cph(MAXVAL(lapw%nv)))
    ALLOCATE(plegend(MAXVAL(lapw%nv),0:atoms%lmaxd))
    ALLOCATE(dplegend(MAXVAL(lapw%nv),0:atoms%lmaxd))
    plegend=0.0
    plegend(:,0)=1.0
    dplegend(:,0)=0.e0
    dplegend(:,1)=1.e0
    !$OMP  DO SCHEDULE(DYNAMIC,1)
    DO  ki =  mpi%n_rank+1, lapw%nv(1), mpi%n_size
       kii=(ki-1)/mpi%n_size+1
       !--->       legendre polynomials
       DO kj = 1,ki
          plegend(kj,1) = DOT_PRODUCT(lapw%gk(:,kj,1),lapw%gk(:,ki,1))
       END DO
       DO l = 1,atoms%lmax(n) - 1
          plegend(:ki,l+1) = fleg1(l)*plegend(:ki,1)*plegend(:ki,l) - fleg2(l)*plegend(:ki,l-1)
          dplegend(:ki,l+1)=REAL(l+1)*plegend(:ki,l)+plegend(:ki,1)*dplegend(:ki,l)
       END DO
       !--->             set up phase factors
       cph = 0.0
       ski = lapw%gvec(:,ki,1)
       DO nn = SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n))
          tnn = tpi_const*atoms%taual(:,nn)
          DO kj = 1,ki
             cph(kj) = cph(kj) +&
                  CMPLX(COS(DOT_PRODUCT(ski-lapw%gvec(:,kj,1),tnn)),&
                  SIN(DOT_PRODUCT(lapw%gvec(:,kj,1)-ski,tnn)))
          END DO
       END DO
       !Set up spinors...
       CALL hsmt_spinor_soc(n,ki,nococonv,lapw,chi,angso)

       !--->          update overlap and l-diagonal hamiltonian matrix
       DO  l = 1,atoms%lmax(n)
          DO j1=1,2
             DO j2=1,2
             !DO j2=j1,j1
                DO kj = 1,ki
                   fct  =cph(kj) * dplegend(kj,l)*fl2p1(l)*(&
                        fj(ki,l,j1)*fj(kj,l,j2) *td%rsoc%rsopp(n,l,j1,j2) + &
                        fj(ki,l,j1)*gj(kj,l,j2) *td%rsoc%rsopdp(n,l,j1,j2) + &
                        gj(ki,l,j1)*fj(kj,l,j2) *td%rsoc%rsoppd(n,l,j1,j2) + &
                        gj(ki,l,j1)*gj(kj,l,j2) *td%rsoc%rsopdpd(n,l,j1,j2)) &
                        * angso(kj,j1,j2)
                   hmat(1,1)%data_c(kj,kii)=hmat(1,1)%data_c(kj,kii) + chi(1,1,j1,j2)*fct
                   hmat(1,2)%data_c(kj,kii)=hmat(1,2)%data_c(kj,kii) + chi(1,2,j1,j2)*fct
                   hmat(2,1)%data_c(kj,kii)=hmat(2,1)%data_c(kj,kii) + chi(2,1,j1,j2)*fct
                   hmat(2,2)%data_c(kj,kii)=hmat(2,2)%data_c(kj,kii) + chi(2,2,j1,j2)*fct
                ENDDO
             ENDDO
          ENDDO
          !--->          end loop over l
       ENDDO
       !--->    end loop over ki
    ENDDO
    !$OMP END DO
    !--->       end loop over atom types (ntype)
    DEALLOCATE(plegend)
    DEALLOCATE(cph)
    !$OMP END PARALLEL
    CALL timestop("offdiagonal soc-setup")

    if (atoms%nlo(n)>0) call hsmt_soc_offdiag_LO(n,atoms,cell,mpi,nococonv,lapw,sym,td,usdus,fj,gj,hmat)

    RETURN
  END SUBROUTINE hsmt_soc_offdiag

  SUBROUTINE hsmt_soc_offdiag_LO(n,atoms,cell,mpi,nococonv,lapw,sym,td,ud,fj,gj,hmat)
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_hsmt_spinor
    USE m_setabc1lo
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_usdus),INTENT(IN)      :: ud
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)!(2,2)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: fj(:,0:,:),gj(:,0:,:)
    !     ..
    !     .. Local Scalars ..
    REAL tnn(3),ski(3)
    INTEGER ki,kj,l,nn,j1,j2,lo,ilo,locol,lorow,ll,nkvec,nkvecp,na,invsfct
    COMPLEX :: fct
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
    COMPLEX:: chi(2,2,2,2),angso(lapw%nv(1),2,2)
    REAL, ALLOCATABLE :: plegend(:,:),dplegend(:,:)
    COMPLEX, ALLOCATABLE :: cph(:)
    REAL                 :: alo1(atoms%nlod,2),blo1(atoms%nlod,2),clo1(atoms%nlod,2)
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

    DO na=sum(atoms%neq(:n-1))+1,sum(atoms%neq(:n))
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
          DO nkvec = 1,invsfct* (2*l+1)
            locol= lapw%nv(1)+lapw%index_lo(lo,na)+nkvec !this is the column of the matrix
            IF (MOD(locol-1,mpi%n_size) == mpi%n_rank) THEN !only this MPI rank calculates this column
              locol=(locol-1)/mpi%n_size+1 !this is the column in local storage
              ki=lapw%kvec(nkvec,lo,na) !this LO is attached to this k+G

              !--->       legendre polynomials
              DO kj = 1,lapw%nv(1)
                plegend(kj,1) = DOT_PRODUCT(lapw%gk(:,kj,1),lapw%gk(:,ki,1))
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
                SIN(DOT_PRODUCT(lapw%gvec(:,kj,1)-ski,tnn)))
              END DO
              !Set up spinors...
              CALL hsmt_spinor_soc(n,ki,nococonv,lapw,chi,angso)

              DO j1=1,2
                DO j2=1,2
                  !DO j2=j1,j1
                  !---> update l-diagonal hamiltonian matrix with LAPW,LO contribution
                  DO kj = 1,lapw%nv(j2)
                    fct  =cph(kj) * dplegend(kj,l)*fl2p1(l)*(&
                    alo1(lo,j1)*fj(kj,l,j2) *td%rsoc%rsopp(n,l,j1,j2) + &
                    alo1(lo,j1)*gj(kj,l,j2) *td%rsoc%rsopdp(n,l,j1,j2) + &
                    blo1(lo,j1)*fj(kj,l,j2) *td%rsoc%rsoppd(n,l,j1,j2) + &
                    blo1(lo,j1)*gj(kj,l,j2) *td%rsoc%rsopdpd(n,l,j1,j2)+ &
                    clo1(lo,j1)*fj(kj,l,j2) *td%rsoc%rsopplo(n,lo,j1,j2) + &
                    clo1(lo,j1)*gj(kj,l,j2) *td%rsoc%rsopdplo(n,lo,j1,j2)) &
                    * angso(kj,j1,j2)
                    hmat(1,1)%data_c(kj,locol)=hmat(1,1)%data_c(kj,locol) + chi(1,1,j1,j2)*fct
                    hmat(1,2)%data_c(kj,locol)=hmat(1,2)%data_c(kj,locol) + chi(1,2,j1,j2)*fct
                    hmat(2,1)%data_c(kj,locol)=hmat(2,1)%data_c(kj,locol) + chi(2,1,j1,j2)*fct
                    hmat(2,2)%data_c(kj,locol)=hmat(2,2)%data_c(kj,locol) + chi(2,2,j1,j2)*fct
                  ENDDO
                  !Update LO-LO part
                  DO ilo=1,atoms%nlo(n)
                    if (l == atoms%llo(ilo,n)) THEN !LO with same L found....
                      DO nkvecp = 1,invsfct* (2*l+1)
                        kj=lapw%kvec(nkvecp,ilo,na) !this LO is attached to this k+G
                        lorow= lapw%nv(1)+lapw%index_lo(ilo,na)+nkvecp !local row
                        fct  =cph(kj) * dplegend(kj,l)*fl2p1(l)*(&
                        alo1(lo,j1)*alo1(ilo,j2) *td%rsoc%rsopp(n,l,j1,j2) + &
                        alo1(lo,j1)*blo1(ilo,j2) *td%rsoc%rsopdp(n,l,j1,j2) + &
                        alo1(lo,j1)*clo1(ilo,j2) *td%rsoc%rsoplop(n,ilo,j1,j2) + &
                        blo1(lo,j1)*alo1(ilo,j2) *td%rsoc%rsoppd(n,l,j1,j2) + &
                        blo1(lo,j1)*blo1(ilo,j2) *td%rsoc%rsopdpd(n,l,j1,j2)+ &
                        blo1(lo,j1)*clo1(ilo,j2) *td%rsoc%rsoplopd(n,ilo,j1,j2)+ &
                        clo1(lo,j1)*alo1(ilo,j2) *td%rsoc%rsopplo(n,lo,j1,j1) + &
                        clo1(lo,j1)*blo1(ilo,j2) *td%rsoc%rsopdplo(n,lo,j1,j1)+ &
                        clo1(lo,j1)*clo1(ilo,j2) *td%rsoc%rsoploplop(n,lo,ilo,j1,j2)) &
                        * angso(kj,j1,j2)
                        hmat(1,1)%data_c(lorow,locol)=hmat(1,1)%data_c(lorow,locol) + chi(1,1,j1,j2)*fct
                        hmat(1,2)%data_c(lorow,locol)=hmat(1,2)%data_c(lorow,locol) + chi(1,2,j1,j2)*fct
                        hmat(2,1)%data_c(lorow,locol)=hmat(2,1)%data_c(lorow,locol) + chi(2,1,j1,j2)*fct
                        hmat(2,2)%data_c(lorow,locol)=hmat(2,2)%data_c(lorow,locol) + chi(2,2,j1,j2)*fct
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDIF !This PE works on LO
          ENDDO!LO
          !--->    end loop over ki
        ENDDO
      ENDIF
    ENDDO
    CALL timestop("offdiagonal soc-setup LO")

    RETURN
  END SUBROUTINE hsmt_soc_offdiag_LO
END MODULE m_hsmt_soc_offdiag
