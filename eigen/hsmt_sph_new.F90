MODULE m_hsmt_sph
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_sph(dimension,atoms,SUB_COMM,n_size,n_rank,sphhar,isp,ab_dim,&
       input,hlpmsize,noco,l_socfirst,cell,nintsp, lapw,el,usdus,vr,gk,rsoc,isigma, aa,bb,fj,gj)

#include"cpp_double.h"
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_sphbes
    USE m_dsphbs
    USE m_ylm
    USE m_hsmt_socinit,ONLY:t_rsoc
    USE m_hsmt_spinor
    USE m_radovlp
#ifdef CPP_MPI
    USE m_mingeselle
#endif
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN):: dimension
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(INOUT)  :: lapw!lpaw%nv_tot is updated
    TYPE(t_usdus),INTENT(INOUT) :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp,ab_dim
    INTEGER, INTENT (IN) :: SUB_COMM,n_size,n_rank 
    INTEGER, INTENT (IN) :: hlpmsize,nintsp
    LOGICAL, INTENT (IN) :: l_socfirst
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntypd,dimension%jspd)
    REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    REAL,INTENT(IN)      :: gk(:,:,:)
    COMPLEX,INTENT(IN)   :: isigma(2,2,3)
    TYPE(t_rsoc),INTENT(IN) :: rsoc


#ifdef CPP_INVERSION
    REAL,    INTENT (INOUT) :: aa(:),bb(:)!(matsize)
#else
    COMPLEX, INTENT (INOUT) :: aa(:),bb(:)
#endif
    REAL,INTENT(OUT) :: fj(:,0:,:,:),gj(:,0:,:,:)
    !     ..
    !     .. Local Scalars ..
    REAL con1,ff,gg,gs,ws,tnn(3), elall,fct,fjkiln,gjkiln,ddnln,ski(3)
    REAL apw_lo1,apw_lo2,apw1,w1

    COMPLEX chi11,chi21,chi22,capw1
    INTEGER ii,iii,ij,k,ki,kj,l,lo,n,n0,n1,nn,kjmax, nsp,iintsp,jintsp
    INTEGER nc ,kii

    !     ..
    !     .. Local Arrays ..
    REAL fb(0:atoms%lmaxd),fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)     
    REAL gb(0:atoms%lmaxd)
   

    REAL, ALLOCATABLE :: plegend(:,:)
    REAL, ALLOCATABLE :: rph(:),cph(:)
  

    !     ..
    con1 = fpi_const/sqrt(cell%omtil)
    DO l = 0,atoms%lmaxd
       fleg1(l) = real(l+l+1)/real(l+1)
       fleg2(l) = real(l)/real(l+1)
       fl2p1(l) = real(l+l+1)/fpi_const
    END DO
    DO iintsp = 1,nintsp
       DO n = 1,atoms%ntype

          DO k = 1,lapw%nv(iintsp)
             gs = lapw%rk(k,iintsp)*atoms%rmt(n)
             CALL sphbes(atoms%lmax(n),gs, fb)
             CALL dsphbs(atoms%lmax(n),gs,fb, gb)
             DO l = 0,atoms%lmax(n)
                !---> set up wronskians for the matching conditions for each ntype
                ws = con1/(usdus%uds(l,n,isp)*usdus%dus(l,n,isp)&
                     - usdus%us(l,n,isp)*usdus%duds(l,n,isp))
                ff = fb(l)
                gg = lapw%rk(k,iintsp)*gb(l)
        !        IF ( apw(l) ) THEN
        !           fj(k,l,n,iintsp) = 1.0*con1 * ff / usdus%us(l,n,isp)
        !           gj(k,l,n,iintsp) = 0.0d0
        !        ELSE
                   !--->                 in a spin-spiral calculation fj and gj are needed
                   !--->                 both interstitial spin directions at the same time
                   !--->                 In all other cases iintsp runs from 1 to 1.
                   fj(k,l,n,iintsp) = ws * ( usdus%uds(l,n,isp)*gg - usdus%duds(l,n,isp)*ff )
                   gj(k,l,n,iintsp) = ws * ( usdus%dus(l,n,isp)*ff - usdus%us(l,n,isp)*gg )
        !        ENDIF
             ENDDO
          ENDDO ! k = 1, lapw%nv
       ENDDO    ! n = 1,atoms%ntype
 
    ENDDO       ! iintsp = 1,nintsp
    !
    ALLOCATE(rph(dimension%nvd))
    ALLOCATE(cph(dimension%nvd))
    ALLOCATE(plegend(dimension%nvd,0:atoms%lmaxd))
  
    plegend=0.0
    plegend(:,0)=1.0
    DO iintsp = 1,nintsp
       DO jintsp = 1,iintsp
          nc = 0                                    ! maybe IF (iintsp.EQ
          lapw%nv_tot = lapw%nv(iintsp)
          DO  kii =  n_rank, lapw%nv_tot-1, n_size
             ki = mod(kii,lapw%nv(iintsp)) + 1
             nc = nc + 1
             !$          nc= 1+ (kii-n_rank)/n_size
             iii = nc*(nc-1)/2 * n_size - (nc-1)*(n_size - n_rank - 1)
             !            ii = nc*(nc+1)/2 * n_size -  nc   *(n_size - n_rank - 1)

             kjmax = ki
             ski = (/lapw%k1(ki,iintsp),lapw%k2(ki,iintsp),lapw%k3(ki,iintsp)/) 
             !--->       legendre polynomials
             DO kj = 1,kjmax
                plegend(kj,1) = dot_product(gk(kj,:,jintsp),gk(ki,:,iintsp))
             END DO
             DO l = 1,maxval(atoms%lmax) - 1
                plegend(:,l+1) = fleg1(l)*plegend(:,1)*plegend(:,l) - fleg2(l)*plegend(:,l-1)
             END DO
             !--->       loop over equivalent atoms
             n1 = 0
             DO n = 1,atoms%ntype

                !---> pk non-collinear
                n0 = n1 + 1
                n1 = n1 + atoms%neq(n)
                rph(:) = 0.0
                cph(:) = 0.0
                DO nn = n0,n1
                   tnn = tpi_const*atoms%taual(:,nn)
                   !--->             set up phase factors
                   DO kj = 1,kjmax
                      rph(kj) = rph(kj) +&
                           cos(dot_product(ski-(/lapw%k1(kj,jintsp),lapw%k2(kj,jintsp),lapw%k3(kj,jintsp)/),tnn))
#ifndef CPP_INVERSION
                      !--->                if the system does not posses inversion symmetry
                      !--->                the complex part of the exponential is needed.
                      cph(kj) = cph(kj) +&
                           sin(dot_product((/lapw%k1(kj,jintsp),lapw%k2(kj,jintsp),lapw%k3(kj,jintsp)/)-ski,tnn))
#endif
                   END DO
                END DO

                !--->          update overlap matrix
                DO kj = 1,ki
                   fct = dot_product(plegend(kj,:atoms%lmax(n))*fl2p1(:atoms%lmax(n)),fj(ki,:atoms%lmax(n),n,iintsp)*fj(kj,:atoms%lmax(n),n,jintsp)&
                        +gj(ki,:atoms%lmax(n),n,iintsp)*gj(kj,:atoms%lmax(n),n,jintsp)*usdus%ddn(:atoms%lmax(n),n,isp))
                   
                   ij = iii + kj
#ifdef CPP_INVERSION
                   bb(ij) = bb(ij) + rph(kj) * fct
#else
                   bb(ij) = bb(ij) + cmplx(rph(kj),cph(kj))*fct
#endif
                END DO
             ENDIF

             
                !--->       end loop over atom types (ntype)
             enddo

             !--->    end loop over ki

          enddo
          !---> end loops over interstitial spins
       ENDDO ! jintsp = 1,iintsp
    ENDDO   ! iintsp = 1,nintsp
    deallocate(plegend)
    deallocate(rph,cph)

    RETURN
  END SUBROUTINE hsmt_sph
END MODULE m_hsmt_sph
