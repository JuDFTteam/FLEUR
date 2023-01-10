!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_offdiag
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_offdiag(n,atoms,fmpi,nococonv,lapw,td,usdus,fjgj,ispin,jspin,iintsp,jintsp,hmat)
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_hsmt_spinor
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_usdus),INTENT(IN)      :: usdus
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_fjgj),INTENT(IN)       :: fjgj
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:)!(2,2)

    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,ispin,jspin,iintsp,jintsp
    !     ..
    !     ..
    !     .. Local Scalars ..
    REAL tnn(3),ski(3)
    INTEGER kii,ki,kj,l,nn,s
    COMPLEX :: fct
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
    REAL fl2p1bt(0:atoms%lmaxd)
    REAL qssbti(3),qssbtj(3)
    COMPLEX:: chi(2,2,2,2)
    REAL, ALLOCATABLE :: plegend(:,:)
    COMPLEX, ALLOCATABLE :: cph(:)

    CALL timestart("offdiagonal setup")

    CALL hsmt_spinor_soc(n,1,nococonv,lapw,chi)



    DO l = 0,atoms%lmaxd
       fleg1(l) = REAL(l+l+1)/REAL(l+1)
       fleg2(l) = REAL(l)/REAL(l+1)
       fl2p1(l) = REAL(l+l+1)/fpi_const
       fl2p1bt(l) = fl2p1(l)*0.5
    END DO
    !$OMP PARALLEL DEFAULT(SHARED)&
    !$OMP PRIVATE(kii,ki,ski,kj,plegend,l)&
    !$OMP PRIVATE(cph,nn,tnn)&
    !$OMP PRIVATE(fct,s)
    ALLOCATE(cph(MAXVAL(lapw%nv)))
    ALLOCATE(plegend(MAXVAL(lapw%nv),0:atoms%lmaxd))
    plegend=0.0
    plegend(:,0)=1.0
    qssbti=MERGE(- nococonv%qss/2,+ nococonv%qss/2,iintsp.EQ.1)
    qssbtj=MERGE(- nococonv%qss/2,+ nococonv%qss/2,jintsp.EQ.1)
    !$OMP  DO SCHEDULE(DYNAMIC,1)
    DO  ki =  fmpi%n_rank+1, lapw%nv(iintsp), fmpi%n_size
       kii=(ki-1)/fmpi%n_size+1
       !--->       legendre polynomials
       DO kj = 1,ki
          plegend(kj,1) = DOT_PRODUCT(lapw%gk(:,kj,jintsp),lapw%gk(:,ki,iintsp))
       END DO
       DO l = 1,atoms%lmax(n) - 1
          plegend(:ki,l+1) = fleg1(l)*plegend(:ki,1)*plegend(:ki,l) - fleg2(l)*plegend(:ki,l-1)
       END DO
       !--->             set up phase factors
       cph = 0.0
       ski = lapw%gvec(:,ki,iintsp) + qssbti
       DO nn = atoms%firstAtom(n), atoms%firstAtom(n) + atoms%neq(n) - 1
          tnn = tpi_const*atoms%taual(:,nn)
          DO kj = 1,ki
             cph(kj) = cph(kj) +&
                  CMPLX(COS(DOT_PRODUCT(ski-lapw%gvec(:,kj,jintsp)+qssbtj,tnn)),&
                  SIN(DOT_PRODUCT(lapw%gvec(:,kj,jintsp)+qssbtj-ski,tnn)))
          END DO
       END DO

       !--->          update overlap and l-diagonal hamiltonian matrix
       s=atoms%lnonsph(n)+1
       DO  l = 0,atoms%lnonsph(n)
          DO kj = 1,ki
             fct  =cph(kj) * plegend(kj,l)*fl2p1(l)*(&
                  fjgj%fj(ki,l,ispin,iintsp)*fjgj%fj(kj,l,jspin,jintsp) *td%h_off(l,l,n,ispin,jspin) + &
                  fjgj%fj(ki,l,ispin,iintsp)*fjgj%gj(kj,l,jspin,jintsp) *td%h_off(l,l+s,n,ispin,jspin) + &
                  fjgj%gj(ki,l,ispin,iintsp)*fjgj%fj(kj,l,jspin,jintsp) *td%h_off(l+s,l,n,ispin,jspin) + &
                  fjgj%gj(ki,l,ispin,iintsp)*fjgj%gj(kj,l,jspin,jintsp) *td%h_off(l+s,l+s,n,ispin,jspin)* sqrt(usdus%ddn(l,n,ispin)*usdus%ddn(l,n,jspin)))
             hmat(1,1)%data_c(kj,kii)=hmat(1,1)%data_c(kj,kii) + CONJG(chi(1,1,iintsp,jintsp)*fct)
             hmat(1,2)%data_c(kj,kii)=hmat(1,2)%data_c(kj,kii) + CONJG(chi(1,2,iintsp,jintsp)*fct)
             hmat(2,1)%data_c(kj,kii)=hmat(2,1)%data_c(kj,kii) + CONJG(chi(2,1,iintsp,jintsp)*fct)
             hmat(2,2)%data_c(kj,kii)=hmat(2,2)%data_c(kj,kii) + CONJG(chi(2,2,iintsp,jintsp)*fct)
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
    CALL timestop("offdiagonal setup")

    RETURN
  END SUBROUTINE hsmt_offdiag
END MODULE m_hsmt_offdiag
