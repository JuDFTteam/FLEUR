!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_offdiag
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_offdiag(n,atoms,mpi,isp,input,noco,cell,iintsp,jintsp,chi,lapw,el,e_shift,usdus,fj,gj,smat,hmat)
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_usdus),INTENT(IN)      :: usdus
    CLASS(t_mat),INTENT(INOUT)     :: smat,hmat
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
    COMPLEX, INTENT(IN)  :: chi
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
    REAL,    INTENT (IN) :: e_shift(input%jspins)
    REAL,    INTENT (IN) :: fj(:,0:,:),gj(:,0:,:)
    !     ..
    !     .. Local Scalars ..
    REAL tnn(3), elall,fct,fjkiln,gjkiln,ddnln,ski(3)
    REAL apw_lo1,apw_lo2,apw1,w1

    COMPLEX capw1
    INTEGER kii,ki,kj,l,nn

    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)     
    REAL fl2p1bt(0:atoms%lmaxd)
    REAL qssbti(3),qssbtj(3)
    REAL, ALLOCATABLE :: plegend(:,:)
    COMPLEX, ALLOCATABLE :: cph(:)
    LOGICAL apw(0:atoms%lmaxd)

    CALL timestart("offdiagonal setup")

    DO l = 0,atoms%lmaxd
       fleg1(l) = REAL(l+l+1)/REAL(l+1)
       fleg2(l) = REAL(l)/REAL(l+1)
       fl2p1(l) = REAL(l+l+1)/fpi_const
       fl2p1bt(l) = fl2p1(l)*0.5
    END DO
    !$OMP PARALLEL DEFAULT(SHARED)&
    !$OMP PRIVATE(kii,ki,ski,kj,plegend,l)&
    !$OMP PRIVATE(cph,nn,tnn,fjkiln,gjkiln)&
    !$OMP PRIVATE(w1,apw_lo1,apw_lo2,ddnln,elall,fct,apw1)&
    !$OMP PRIVATE(capw1) 
    ALLOCATE(cph(MAXVAL(lapw%nv)))
    ALLOCATE(plegend(MAXVAL(lapw%nv),0:atoms%lmaxd))
    plegend=0.0
    plegend(:,0)=1.0
    qssbti=MERGE(- noco%qss/2,+ noco%qss/2,iintsp.EQ.1)
    qssbtj=MERGE(- noco%qss/2,+ noco%qss/2,jintsp.EQ.1)
    !$OMP  DO SCHEDULE(DYNAMIC,1)
    DO  ki =  mpi%n_rank+1, lapw%nv(iintsp), mpi%n_size
       kii=(ki-1)/mpi%n_size+1
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
       DO nn = SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n))
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
          ddnln =  usdus%ddn(l,n,isp)
          DO kj = 1,ki
             fct  =cph(kj) * plegend(kj,l)*fl2p1(l)*(&
                  fj(ki,l,iintsp)*fj(kj,l,jintsp) *td%h_off(l,l,n,isp) + &
                  fj(ki,l,iintsp)*gj(kj,l,jintsp) *td%h_off(l,l+s,n,isp) + &
                  gj(ki,l,iintsp)*fj(kj,l,jintsp) *td%h_off(l+s,l,n,isp) + &
                  gj(ki,l,iintsp)*gj(kj,l,jintsp) *td%h_off(l+s,l+s,n,isp)*ddnln)
             hmat(1,1)%data_c(kj,kii)=hmat(1,1)%data_r(kj,kii) + chi(1,1)*fct 
             hmat(1,2)%data_c(kj,kii)=hmat(1,2)%data_r(kj,kii) + chi(1,2)*fct 
             hmat(2,1)%data_c(kj,kii)=hmat(2,1)%data_r(kj,kii) + chi(2,1)*fct 
             hmat(2,2)%data_c(kj,kii)=hmat(2,2)%data_r(kj,kii) + chi(2,2)*fct 
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
